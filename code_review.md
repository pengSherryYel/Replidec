# RepliDec code review checklist

Review the complete pull request diff and verify claims in the pull request
description against the code and check results. Prioritize correctness,
scientific validity, and regression risk over formatting preferences.

Codex should post only concrete merge-blocking P0 or P1 findings. Human
reviewers may also leave clearly labelled, non-blocking P2 improvements or P3
suggestions.

## Severity

- **P0 — critical:** security or credential exposure, data loss, invalid release
  control, or a defect that can silently invalidate scientific results.
- **P1 — high:** a concrete correctness, regression, compatibility, packaging,
  or reproducibility problem that should block merging.
- **P2 — medium:** a worthwhile non-blocking maintainability or test improvement.
- **P3 — low:** an optional suggestion or style preference.

## 1. Architecture and scope

- [ ] The change has one clear purpose and contains no unrelated refactoring.
- [ ] The change fits the existing architecture, or an intentional design
      change and its migration impact are explained explicitly.
- [ ] Public CLI behavior, file formats, defaults, and output locations remain
      backward compatible unless the breaking change is documented.

## 2. Correctness and scientific behavior

- [ ] Inputs, empty data, malformed records, missing files, and external-tool
      failures cannot silently produce misleading results.
- [ ] Biological classification logic, thresholds, and database assumptions
      change only when the pull request explicitly intends and validates them.
- [ ] Results remain deterministic and reproducible for the same input,
      configuration, dependency set, and database version.

## 3. Regression protection

- [ ] New or changed behavior has an automated test at the appropriate level.
- [ ] Affected single-sample and multi-sample workflows remain functional.
- [ ] Error paths and boundary cases touched by the change are tested.
- [ ] Required Flake8 and integration checks pass on the latest commit.

## 4. Packaging and release safety

- [ ] Package version, dependencies, package data, entry points, and supported
      Python versions remain internally consistent.
- [ ] Production automation can act only on a reviewed commit from `main`
      through the protected version-tag workflow.
- [ ] No credentials, tokens, private URLs, unpublished datasets, generated
      files, or local artifacts are committed.

## 5. Review output

For each finding, include:

1. severity and a concise title;
2. the relevant file and smallest useful line range;
3. the concrete failure or regression scenario; and
4. a suggested correction when one is clear.

Do not report speculative or style-only comments as defects. If there are no
actionable findings, say so explicitly and identify any material testing
limitation. Automated review never substitutes for the repository owner's
final human approval.
