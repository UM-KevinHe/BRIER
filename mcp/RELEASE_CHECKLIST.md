# Release checklist - pushing BRIER MCP v1.0.0 to GitHub

This checklist covers everything I (Claude) can't do for you: the actual `git` and GitHub steps. Run from `~/Desktop/Dissertation/BRIER_software/BRIER_MCP/` on your Mac.

---

## 0. Confirm what you're shipping

Run the test suite once more to be sure the working tree is green:

```bash
cd mcp
uv run tests/test_v133.py     # latest specific suite
# Optionally all of them:
for t in tests/test_*.py; do echo "=== $t ==="; uv run "$t"; done
cd ..
```

You want every suite to end in `ALL TESTS PASSED`.

Sanity-check the version:

```bash
python3 -c "import json; print(json.load(open('mcp/manifest.json'))['version'])"
# should print: 1.0.0
```

---

## 1. Decide: subdirectory of UM-KevinHe/BRIER, or its own repo

You picked "subdirectory of UM-KevinHe/BRIER" in our build conversation. That means the MCP server ends up at `UM-KevinHe/BRIER/mcp/` on GitHub. The instructions below assume that path.

If you change your mind and want a separate repo (e.g. `UM-KevinHe/BRIER-MCP`), the steps are similar but you'd `git init` at `mcp/` itself rather than the parent.

---

## 2. Get a working copy of the UM-KevinHe/BRIER repo locally

If you don't already have a clone:

```bash
cd ~/Desktop/Dissertation
git clone https://github.com/UM-KevinHe/BRIER.git BRIER-upstream
cd BRIER-upstream
```

If you already maintain a clone, just `cd` into it.

Verify you have push access:

```bash
git remote -v
# origin should be your fork or the canonical repo
```

If you don't have push access, you'll need to fork on GitHub and add your fork as a remote, or coordinate with whoever owns `UM-KevinHe`.

---

## 3. Copy the v1.0.0 `mcp/` into the upstream tree

If the upstream repo doesn't yet have an `mcp/` folder:

```bash
cp -r ~/Desktop/Dissertation/BRIER_software/BRIER_MCP/mcp ./mcp
```

If it already has one (e.g. from an earlier preview), replace it:

```bash
# from inside the BRIER-upstream working tree
rm -rf mcp
cp -r ~/Desktop/Dissertation/BRIER_software/BRIER_MCP/mcp ./mcp
```

Then make sure no stray cache files snuck in:

```bash
find mcp -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null
find mcp -name ".pytest_cache" -type d -exec rm -rf {} + 2>/dev/null
find mcp -name "*.pyc" -delete 2>/dev/null
```

### 3a. Add `mcp/` to `.Rbuildignore` (important)

The repo root is an R package. When someone runs `devtools::install_github("UM-KevinHe/BRIER")` or `R CMD build .`, R packages everything in the repo root by default. You need to tell R to skip `mcp/`, or:

- The R package tarball will be larger than it should be
- `R CMD check` will warn about non-standard files at the top level
- If the package is ever submitted to CRAN, it would be rejected for unexpected files

The fix is one line:

```bash
# from inside the BRIER repo root
cat .Rbuildignore 2>/dev/null   # see what's there; create if missing
echo '^mcp$' >> .Rbuildignore
```

This is a regex meaning "the folder named exactly `mcp` at the repo root." After this, `mcp/` is invisible to R's build tooling. The folder still lives in the GitHub repo and on every clone; it's just excluded when R bundles the package.

You can verify it works:

```bash
R CMD build . 2>&1 | grep -i "mcp"   # should print nothing
```

(That requires R installed locally and may be slow; the `echo` line alone is enough if you trust the regex.)

Check what'll be committed:

```bash
git status
git diff --stat | tail
```

You should see `mcp/` (and `.Rbuildignore`) in the list. Nothing else should be touched (in particular, the R package source under `R/`, `src/`, `man/`, `DESCRIPTION` should be unchanged).

---

## 4. Commit and push to a branch

Don't push directly to `main` for a first release - open a PR so the repo owners can review.

```bash
git checkout -b mcp/v1.0.0
git add mcp/ .Rbuildignore
git commit -m "Add MCP server v1.0.0 for Claude Desktop integration

First public-stable release of the BRIER MCP. 27 tools, 21 test suites
passing. Lets Claude Desktop drive the BRIER R package via natural
language: inspect data, recommend the right fit family, run fitting and
selection, generate HTML reports. All computation runs locally through
the user's R installation; data never leaves the machine.

See mcp/README.md for install instructions and mcp/RELEASE_NOTES_v1.0.0.md
for the full release writeup."

git push origin mcp/v1.0.0
```

Then open a PR on GitHub from `mcp/v1.0.0` into `main`.

---

## 5. After the PR is merged, tag the release

Once merged to `main`:

```bash
git checkout main
git pull origin main

# Create an annotated tag.
# We use the prefix `mcp-` so MCP releases don't collide with future
# tags for the BRIER R package itself (which would be `v0.3.0` etc.
# unprefixed). Mirror of the BregSurv convention (their MCP tag is
# `mcpb-v1.0.0`).
git tag -a mcp-v1.0.0 -m "BRIER MCP v1.0.0 - first public-stable release"
git push origin mcp-v1.0.0
```

---

## 6. Create the GitHub Release

On the repo page → **Releases** → **Draft a new release**.

- **Tag**: `mcp-v1.0.0` (the one you just pushed)
- **Title**: `BRIER MCP v1.0.0`
- **Description**: copy-paste the contents of `mcp/RELEASE_NOTES_v1.0.0.md`
- **Assets**: optional - upload a source tarball if you want non-git users to be able to download a zip:
  ```bash
  cd /path/to/BRIER-upstream
  git archive --format=tar.gz --prefix=BRIER-mcp-1.0.0/ -o ~/Desktop/BRIER-mcp-1.0.0.tar.gz mcp-v1.0.0 mcp/
  ```
  Then drag the resulting `BRIER-mcp-1.0.0.tar.gz` into the GitHub release form's assets area.
- Check **"Set as the latest release"** if this is the project's first MCP release; otherwise judgment call.

Click **Publish release**.

---

## 7. (Optional) Announce somewhere

If the BRIER project has a website (`um-kevinhe.github.io/BRIER/`), consider adding an MCP install section that links to:

- `https://github.com/UM-KevinHe/BRIER/blob/main/mcp/README.md` (install instructions)
- `https://github.com/UM-KevinHe/BRIER/releases/tag/mcp-v1.0.0` (the release page)

---

## 8. For future releases (v1.1.0+)

The pattern is the same, just shorter:

1. Develop on a branch (e.g. `mcp/v1.1.0`).
2. Update `mcp/manifest.json` version + `mcp/pyproject.toml` version.
3. Append a new section to `mcp/DEVLOG.md` describing what changed.
4. Update the **Status** block at the top of `mcp/README.md`.
5. Create a `mcp/RELEASE_NOTES_v1.X.0.md` (copy v1.0.0 as a template).
6. PR → merge → tag with `mcp-` prefix (e.g. `mcp-v1.1.0`) → GitHub Release.

**Always prefix MCP tags with `mcp-`** so they sort separately from R-package tags. The R package presumably uses bare tags (`v0.3.0`); MCP uses `mcp-v1.X.Y`.

If breaking changes (renaming tools, removing args, changing env vars), bump to v2.0.0 instead of v1.X.0 and call it out clearly in the release notes.

---

## Things I (Claude) can't do for you

- `git push` (your machine, your credentials).
- Create the GitHub PR / Release.
- Confirm your push access on `UM-KevinHe/BRIER`.
- Verify the install path actually works on your Mac (you have to do this by following `mcp/README.md` yourself once and reporting any issues).

If anything in this checklist is unclear or breaks, paste the error back and I'll help debug.
