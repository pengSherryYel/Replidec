# Contributing to RepliDec

RepliDec uses a deliberately small two-branch workflow. Development should stay
easy for contributors, while every production change remains reviewable.

## Branch roles

| Branch | Purpose | Working rule |
| --- | --- | --- |
| `develop` | Shared integration branch for ongoing work | Contributors may push directly. Feature branches are optional for larger or experimental changes. |
| `main` | Reviewed, releasable history | Never push directly. Changes enter only through a reviewed pull request from `develop`. |

Formal version tags and production package releases are created from `main`
only. There is no fixed weekly release day: open a `develop` to `main` pull
request when the accumulated changes form a coherent release candidate.

## Activating the workflow

The repository owner initializes `develop` once, after the CI/release guardrails
and this workflow documentation have both been reviewed and merged into `main`:

```bash
git switch main
git pull --ff-only origin main
git switch --create develop
git push --set-upstream origin develop
```

Create `develop` from the current reviewed `main`, not from an unmerged pull
request branch. Until `origin/develop` exists, contributors should continue to
use focused pull-request branches targeting `main`; the day-to-day workflow
below becomes active after the owner completes this one-time initialization.

## Day-to-day development

1. Start from the current shared branch:

   ```bash
   git fetch origin
   git switch develop
   git pull --ff-only origin develop
   ```

2. For a focused change, work directly on `develop`. For a risky, large, or
   collaborative change, create a feature branch and merge it back into
   `develop` when ready.
3. Keep `develop` usable. Do not force-push it or rewrite its published
   history, and pull before pushing to reduce avoidable conflicts.
4. Run the relevant local checks before pushing.

## Local lint setup

Install the development tools once:

```bash
python -m pip install --requirement requirements-dev.txt
pre-commit install --hook-type pre-push
```

The pre-push hook runs Flake8 automatically. It can also be run explicitly:

```bash
python -m flake8 .
```

Run tests relevant to the changed behavior. The full integration workflow,
including the reference database and external tools, runs in GitHub Actions.

## Promoting `develop` to `main`

Create a pull request with `develop` as the source and `main` as the target.
Complete the pull request template with:

- the purpose and scope of the change;
- any architecture, scientific-method, CLI, output, or compatibility impact;
- the exact validation performed;
- concrete regression risks and suggested review focus;
- whether the change should produce a formal release.

The pull request may merge only after:

1. required lint and integration checks pass on the latest commit;
2. Codex automated review findings are resolved or answered;
3. review conversations are resolved; and
4. the repository owner completes the final human review and approval.

New commits after approval require the checks and human review to be repeated.
AI-assisted review is supporting evidence, not permission to merge.

After the pull request merges, synchronize `develop` with `main` before the
next development cycle. The repository owner separately decides whether and
when to create a protected `vX.Y.Z` release tag.
