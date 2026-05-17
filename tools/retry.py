"""Exponential-backoff retry wrapper for flaky network calls (MSA server, HF Hub).

Brief-locked retry schedule: 60 s, 300 s, 900 s (3 attempts after the first).
We do NOT add jitter, since the brief mandates exact delays.

Per-call timeout is enforced via a SIGALRM-based watchdog (POSIX only).
"""
from __future__ import annotations

import logging
import signal
import time
from dataclasses import dataclass
from typing import Callable, Iterable, TypeVar

T = TypeVar("T")

log = logging.getLogger("tools.retry")

# Brief-locked schedule. First attempt has no preceding delay.
DEFAULT_BACKOFF_SECONDS: tuple[int, ...] = (60, 300, 900)
DEFAULT_PER_CALL_TIMEOUT_SECONDS: int = 2 * 60 * 60  # 2h per prediction


class CallTimeoutError(TimeoutError):
    """Raised when a single attempt exceeds `timeout_seconds`."""


@dataclass
class RetryStats:
    """Bookkeeping for one retry session — useful for logging."""
    attempts: int = 0
    elapsed_seconds: float = 0.0
    exceptions: list[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.exceptions is None:
            self.exceptions = []


def _alarm_handler(signum, frame):
    raise CallTimeoutError(f"call exceeded per-attempt timeout")


def call_with_timeout(fn: Callable[[], T], timeout_seconds: int) -> T:
    """Run `fn()` and raise CallTimeoutError after `timeout_seconds`.

    Uses SIGALRM, which is process-wide. Safe to nest only if outer caller does
    not also install a SIGALRM handler. For our use case (ColabFold MSA → AF2
    prediction), this is the only timer in flight.
    """
    if timeout_seconds <= 0:
        return fn()
    prev_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(timeout_seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev_handler)


def retry(
    fn: Callable[[], T],
    *,
    backoff_seconds: Iterable[int] = DEFAULT_BACKOFF_SECONDS,
    per_call_timeout_seconds: int = DEFAULT_PER_CALL_TIMEOUT_SECONDS,
    retry_on: tuple[type[BaseException], ...] = (Exception,),
    label: str = "call",
    sleep_fn: Callable[[float], None] = time.sleep,
    stats: RetryStats | None = None,
) -> T:
    """Run `fn` with exponential backoff retries.

    Sequence is: attempt 1 → on failure wait backoff_seconds[0], attempt 2 → wait
    backoff_seconds[1], attempt 3 → wait backoff_seconds[2], attempt 4 → propagate.
    With the default schedule (60, 300, 900), the maximum delay before giving up
    is ~21.3 minutes plus per-attempt timeouts.

    Each individual attempt is wrapped in `call_with_timeout` and will raise
    CallTimeoutError if it does not finish within `per_call_timeout_seconds`.
    CallTimeoutError is retried like any other Exception by default.

    Args:
        fn:                          zero-arg callable to run.
        backoff_seconds:             post-failure delays. With N values, total
                                     attempts = N + 1.
        per_call_timeout_seconds:    SIGALRM-based timeout for each attempt.
                                     Pass 0 or negative to disable.
        retry_on:                    exception classes that trigger a retry.
                                     Anything else propagates immediately.
        label:                       short string for log lines.
        sleep_fn:                    inject for testing.
        stats:                       optional RetryStats to populate.
    """
    backoffs = list(backoff_seconds)
    if stats is None:
        stats = RetryStats()
    start = time.monotonic()
    last_exc: BaseException | None = None
    for attempt_idx in range(len(backoffs) + 1):
        stats.attempts = attempt_idx + 1
        try:
            result = call_with_timeout(fn, per_call_timeout_seconds)
            stats.elapsed_seconds = time.monotonic() - start
            log.info("%s: ok on attempt %d", label, stats.attempts)
            return result
        except retry_on as e:  # noqa: PERF203 — we want the explicit catch
            last_exc = e
            stats.exceptions.append(f"{type(e).__name__}: {e}")
            log.warning("%s: attempt %d failed (%s: %s)",
                        label, stats.attempts, type(e).__name__, e)
            if attempt_idx >= len(backoffs):
                break
            delay = backoffs[attempt_idx]
            log.info("%s: sleeping %ds before attempt %d",
                     label, delay, stats.attempts + 1)
            sleep_fn(delay)
    stats.elapsed_seconds = time.monotonic() - start
    assert last_exc is not None
    raise last_exc
