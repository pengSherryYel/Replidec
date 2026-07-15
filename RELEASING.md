# RepliDec release process

RepliDec uses GitHub Actions to test every pull request and release candidate.
Production releases publish to PyPI, create a GitHub Release, and open an
update pull request in `bioconda/bioconda-recipes`. Bioconda publishes both the
Conda package and its BioContainer after that pull request passes Bioconda CI
and is merged.

## One-time repository setup

1. Enable branch protection or a ruleset for `main`:
   - require pull requests;
   - require review from the code owner;
   - dismiss stale approvals when new commits are pushed;
   - require the `Flake8` status check;
   - require the `RepliDec integration` status check;
   - require branches to be up to date before merging;
   - disallow force pushes and branch deletion;
   - enable GitHub auto-merge if desired.
2. Protect release tags matching `v*` from deletion and force updates.
3. Create a GitHub environment named `pypi`.
4. In the PyPI Replidec project, configure a Trusted Publisher for:
   - owner: `pengSherryYel`;
   - repository: `Replidec`;
   - workflow: `release.yml`;
   - environment: `pypi`.
5. Fork `bioconda/bioconda-recipes` to the account that will submit recipe
   updates. Create a GitHub environment named `bioconda` and add:
   - secret `BIOCONDA_PR_TOKEN`: a token able to push to that fork and open a
     pull request against `bioconda/bioconda-recipes`;
   - optional variable `BIOCONDA_FORK_OWNER` when the fork owner differs from
     the RepliDec repository owner.

No Anaconda token or Docker Hub credential is required. Bioconda owns the
final Conda upload and creates the versioned BioContainer on Quay.

## Prepare a release

1. Update `__version__` in `Replidec/__init__.py`. This is the single source
   used by packaging and the CLI.
2. Open a pull request and wait for `RepliDec integration` to pass.
3. Optionally run the `Release` workflow manually with `ref` set to the
   candidate branch and `tag` set to the expected `vX.Y.Z` value. Manual runs
   are always dry runs and cannot change an external service.
4. Merge the release pull request.
5. Create and push an annotated tag that exactly matches the package version:

   ```bash
   git tag -a v0.3.6 -m "RepliDec 0.3.6"
   git push origin v0.3.6
   ```

The tag starts the production workflow. PyPI versions are immutable, so never
reuse or force-update a release tag.

## After the workflow

The workflow summary contains the Bioconda pull request URL. Bioconda review
and merge remain external quality gates. After the recipe is merged, users can
install the new Conda package and pull the automatically generated container:

```bash
conda install -c conda-forge -c bioconda replidec
docker pull quay.io/biocontainers/replidec:<bioconda-tag>
```
