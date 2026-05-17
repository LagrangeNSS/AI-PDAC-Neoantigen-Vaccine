"""Unit tests for tools.retry."""
from __future__ import annotations

import time

import pytest

from tools.retry import (
    CallTimeoutError,
    DEFAULT_BACKOFF_SECONDS,
    RetryStats,
    call_with_timeout,
    retry,
)


def test_default_schedule_is_brief_locked() -> None:
    assert DEFAULT_BACKOFF_SECONDS == (60, 300, 900)


def test_succeeds_on_first_attempt() -> None:
    calls = []

    def fn() -> int:
        calls.append(1)
        return 42

    stats = RetryStats()
    out = retry(fn, sleep_fn=lambda _: None, stats=stats)
    assert out == 42
    assert stats.attempts == 1
    assert calls == [1]


def test_succeeds_after_two_failures() -> None:
    calls = []
    sleeps: list[float] = []

    def fn() -> str:
        calls.append("call")
        if len(calls) < 3:
            raise ConnectionError(f"fail {len(calls)}")
        return "ok"

    stats = RetryStats()
    out = retry(fn, sleep_fn=sleeps.append, stats=stats)
    assert out == "ok"
    assert stats.attempts == 3
    assert sleeps == [60, 300]
    assert len(stats.exceptions) == 2


def test_gives_up_after_full_schedule_and_reraises_last() -> None:
    sleeps: list[float] = []
    counter = {"n": 0}

    def fn() -> None:
        counter["n"] += 1
        raise ValueError(f"attempt {counter['n']}")

    stats = RetryStats()
    with pytest.raises(ValueError, match="attempt 4"):
        retry(fn, sleep_fn=sleeps.append, stats=stats)
    assert stats.attempts == 4
    assert sleeps == [60, 300, 900]


def test_non_retryable_exception_propagates_immediately() -> None:
    calls = []

    def fn() -> None:
        calls.append(1)
        raise KeyboardInterrupt()  # not in retry_on

    with pytest.raises(KeyboardInterrupt):
        retry(fn, sleep_fn=lambda _: None, retry_on=(ConnectionError,))
    assert calls == [1]


def test_call_with_timeout_completes_normally() -> None:
    assert call_with_timeout(lambda: 7, timeout_seconds=2) == 7


def test_call_with_timeout_fires_on_overrun() -> None:
    def slow() -> None:
        time.sleep(3)

    with pytest.raises(CallTimeoutError):
        call_with_timeout(slow, timeout_seconds=1)


def test_per_call_timeout_is_retryable_by_default() -> None:
    calls = []

    def slow_then_fast() -> str:
        calls.append(1)
        if len(calls) == 1:
            time.sleep(3)  # blows the 1 s timeout below
        return "ok"

    out = retry(
        slow_then_fast,
        backoff_seconds=(0,),
        per_call_timeout_seconds=1,
        sleep_fn=lambda _: None,
    )
    assert out == "ok"
    assert calls == [1, 1]


def test_custom_short_schedule() -> None:
    sleeps: list[float] = []
    n = {"i": 0}

    def fn() -> int:
        n["i"] += 1
        if n["i"] < 2:
            raise ConnectionError("boom")
        return 1

    out = retry(fn, backoff_seconds=(5,), sleep_fn=sleeps.append)
    assert out == 1
    assert sleeps == [5]
