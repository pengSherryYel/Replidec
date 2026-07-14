# RepliDec repository guidance

## Working rules

- Treat `develop` as the shared development branch and `main` as protected,
  reviewed, and releasable history.
- Never merge into `main`, create a formal version tag, or start a production
  release without explicit approval from the repository owner.
- Keep changes focused and preserve existing behavior unless the task clearly
  calls for an intentional, documented change.

## Validation

- Before proposing a push or merge, run `python -m flake8 .`.
- Run tests that exercise the changed behavior and report the exact commands
  and results. Use the GitHub integration check for the full external-tool and
  reference-database workflow.
- Do not weaken, skip, or remove a required check merely to obtain a passing
  result.

## Review guidelines

For every pull request review, read and follow `code_review.md`.

Use the pull request description as context, but verify its claims against the
complete diff and the latest check results. Focus automated review comments on
concrete P0 and P1 defects. Codex review is an additional quality check, not
final approval; merging into `main` still requires human confirmation from the
repository owner.
