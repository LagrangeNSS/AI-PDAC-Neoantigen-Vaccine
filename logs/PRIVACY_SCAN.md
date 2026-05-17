# Privacy scan log

Scan performed during the open-source bundling step, before `git init`. The
goal is to make sure no API keys, OAuth tokens, personal credentials,
usernames, or local home/Windows paths leak into the public repository.

## Scans performed

All scans target the staged tree under `./GithubPublish/` only. No file
outside this tree is modified.

### 1. API keys and tokens

Pattern (case-sensitive, ripgrep extended regex):

```
(sk-ant-|sk-proj-|ANTHROPIC_API_KEY|Bearer\s+[A-Za-z0-9]|ghp_|gho_|github_pat_|api[_-]?key\s*=)
```

Result: **0 matches.**

### 2. Generic secret-looking key/value pairs

Pattern (case-insensitive):

```
(token|password|secret|credential)\s*[=:]
```

Result: **0 matches.**

### 3. Usernames and home / Windows paths

The host username detected via `whoami` is `NSS`. Patterns searched:

```
NSS@
/home/NSS/
/mnt/c/Users/NSS/
C:\Users\NSS
/home/<any-username>/
/Users/<any-username>/
the literal string "NSS"
```

Result: **0 matches.** The project root is `/mnt/d/AssignmentOPUS/`,
which is a drive-letter-only path with no user-identifying segment, so
the host username does not propagate into the logs, code, or
documentation.

## Redactions applied

None. The scan returned no hits, so no file content was modified.

## Files scanned

The scans walked every file under `./GithubPublish/`:

- `prompts/RESEARCH_PROMPT.md`
- `tools/**` (Python modules + tests)
- `data/sequences/*.fasta`
- `data/structures/**/*.pdb`
- `results/*.json`
- `figures/*.png` + `figures/*.txt`
- `logs/*.md` (whitelisted set)
- `refs/citations.bib`, `refs/notes.md`
- `RESEARCH_README.md`

## Conclusion

The `./GithubPublish/` tree is clear of secrets and personal identifiers
as far as this scan can determine. Note that this is a regex-based
sanity check, not a guarantee — the repository owner should still
visually skim the bundle once before pushing.
