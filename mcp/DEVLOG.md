# BRIER MCP server (developer notes)

A FastMCP server that exposes BRIER R functions as MCP tools so Claude
Desktop (or any MCP client) can fit transfer-learning genetic risk
prediction models on the user's machine without round-tripping through
hand-written R code.

> **Status: v0.13.3.** Single-default change. `bootstrap_n` in
> `summarize_fit` and `brier_plot_box`, and `replications` in
> `brier_plot_importance`, drop from 100 to 20. Reason: on real
> 10k-SNP fits each bootstrap replicate refits the model, and 100
> replicates ran tens of minutes per plot. n=20 is enough for a quick
> exploratory comparison; users wanting publication-quality variance
> estimates pass `bootstrap_n=100` (or higher) explicitly. The
> `brier_plot_eta` `bootstrap_n` stays at 100 because that plot only
> bootstraps when `bootstrap=True` is set explicitly (deliberate
> opt-in, not on the hot path). 27 tools, 21 test suites passing.

## v0.13.2 (unchanged from this release)

> **v0.13.2** added two wizard-text refinements: the welcome-URL
> append anti-pattern called out explicitly in `_display_instructions`,
> and an EXPRESSION VALIDATOR section in `ai_instructions` correcting
> the assistant's repeated claim that `::` is blocked (it isn't; the
> whitelist of `BRIER::` / `base::` / `stats::` / `utils::` /
> `Matrix::` has been in place since v0.7).

## v0.13.1 (unchanged from this release)

> **v0.13.1** introduced soft bootstrap-plot degradation (default
> `plot_timeout_seconds=None`, fallback snippet only for retry-worthy
> failures, legacy `bootstrap_plot_max_p` kept for back-compat) and
> the ETA GRID DO-NOT-HAND-WRITE wizard guidance.

## v0.13 (unchanged from this release)

> **v0.13** introduced `prep_data` (9 composable operations, audit log,
> summarize_fit integration) and `prep_data_log`.

## v0.12 (unchanged from this release)

> **v0.12** added strict-equality escalation trigger and single-shot
> de-escalation to `brier_auto_tune_eta`. Default
> `near_boundary_top_fraction=0.0`; default `de_escalation_threshold=1.0`,
> `de_escalation_ceiling=2.0`. Old top-20% behavior available by
> passing `near_boundary_top_fraction=0.20`.

## v0.11 (unchanged from this release)

> **v0.11** introduced `brier_auto_tune_eta`: automatic eta-ceiling
> escalation across a fixed ladder (default `[30, 50, 100]`). Plain
> selection tools remain diagnostic-only; auto-escalation is opt-in
> by tool name. Returns full provenance (`final_fit_id`,
> `final_selection_id`, `escalation_history`).

## v0.10.3 features (unchanged)

> **v0.10.3** brought four features: (1) **Principled eta.list default**:
> when `eta_list` is not supplied, the fit tools build a log-spaced grid
> `c(0, exp(seq(log(eta_floor), log(eta_ceiling), length.out=eta_n)))`
> with defaults `eta_floor=0.1`, `eta_ceiling=10`, `eta_n=10`. The
> explicit grid is stored in the fit cache and emitted in reproduce.R
> for exact reproduction. (2) **Boundary-optimum diagnostic**: after
> selection, if the chosen eta sits at the top of the grid, a
> `_notice_eta_boundary` is attached suggesting a refit with a higher
> ceiling. (3) **`brier_plot_selection` tool**: plots the selection
> criterion vs eta directly from the cached selection (no test data
> required). (4) **Per-call `output_dir` override** on `summarize_fit`,
> all three plot tools, and `brier_predict`.

## Directory layout

```
mcp/
├── server.py                    FastMCP server: @mcp.tool definitions, Rscript bridge, wizard.
├── manifest.json                MCPB extension manifest (manifest_version=0.4).
├── pyproject.toml               Runtime deps for the bundle (uv reads this at install).
├── .mcpbignore                  Excludes tests/cache from `mcpb pack`.
├── r_scripts/
│   ├── _common.R                Shared utilities sourced by every dispatcher.
│   ├── inspect_data.R           Describe an .rda/.RData/.rds file's structure.
│   ├── inspect_user_data.R      Heuristic inspection with format support.
│   ├── list_data_directory.R    List .rda/.RData/.rds files in a directory.
│   ├── brier_i.R                Fit BRIERi() (pretrained external + individual target).
│   ├── brier_i_cv.R             Cross-validation tuning for BRIERi.
│   ├── brier_i_selection.R      IC- or validation-set selection on a brier_i fit.
│   ├── brier_full.R             Fit BRIERfull() (pooled-cohort, raw external).
│   ├── brier_full_selection.R   Validation-set selection on a brier_full fit.
│   ├── brier_s.R                Fit BRIERs() (summary-statistics target).
│   ├── brier_s_selection.R      IC- or validation-set selection on a brier_s fit.
│   ├── get_ldb.R                Return Berisa-Pickrell LD block coordinates.
│   ├── cal_ld.R                 Build an LD matrix from a reference panel.
│   ├── brier_predict.R          Predict from any cached fit / selection on new X.
│   └── brier_evaluate.R         Score any cached fit / selection on (X, y).
├── tests/
│   ├── test_inspect_data.py     v0.1 baseline.
│   ├── test_v02.py              v0.2 BRIERi family.
│   ├── test_v03.py              v0.3 predict + evaluate loop.
│   ├── test_v04.py              v0.4 BRIERfull family.
│   ├── test_v05.py              v0.5 BRIERs family + LD utilities.
│   └── test_v06.py              v0.6 wizard.
└── README.md                    This file.
```

Note: `start_analysis` (the wizard) is implemented entirely in
`server.py` as a pure-Python tool. It returns a structured dict with
framing, primers, problem-description questions, routing logic with
size-based recommendations, model paths, preprocessing hints, family
caveats, and reproducibility advice. No R dispatcher needed.

## Local dev workflow (without Claude Desktop)

Fastest iteration loop while building out the tool surface:

1. Install R + the BRIER R package + jsonlite.

   ```r
   install.packages("jsonlite")
   if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
   remotes::install_github("UM-KevinHe/BRIER")
   ```

2. Install `uv` (https://github.com/astral-sh/uv), which bundles Python deps.

3. From this directory, run all twelve smoke tests:

   ```bash
   uv run tests/test_inspect_data.py    # v0.1 baseline (6 PASS)
   uv run tests/test_v02.py             # v0.2 BRIERi family (~19 PASS)
   uv run tests/test_v03.py             # v0.3 predict + evaluate (~22 PASS)
   uv run tests/test_v04.py             # v0.4 BRIERfull on Data_BRIERfull (~22 PASS)
   uv run tests/test_v05.py             # v0.5 BRIERs + LD utilities on Data_BRIERs
   uv run tests/test_v06.py             # v0.6 wizard structure (~60 PASS)
   uv run tests/test_v07.py             # v0.7 data-first flow (~40 PASS)
   uv run tests/test_v071.py            # v0.7.1 fixes (~25 PASS)
   uv run tests/test_v080.py            # v0.8.0 refinements (~30 PASS)
   uv run tests/test_v081.py            # v0.8.1 reverts + preprocess (~25 PASS)
   uv run tests/test_v090.py            # v0.9 plot wrappers (~25 PASS)
   uv run tests/test_v100.py            # v0.10 summarize_fit report (~30 PASS)
   uv run tests/test_v101.py            # v0.10.1 wizard UX cleanup (~25 PASS)
   uv run tests/test_v102.py            # v0.10.2 multi-file data_paths (~25 PASS)
   uv run tests/test_v103.py            # v0.10.3 eta default + selection plot + output_dir + boundary notice (~30 PASS)
   uv run tests/test_v110.py            # v0.11 brier_auto_tune_eta escalation (~25 PASS)
   uv run tests/test_v120.py            # v0.12 strict trigger + de-escalation (~25 PASS)
   uv run tests/test_v130.py            # v0.13 prep_data + audit log + integration (~40 PASS)
   uv run tests/test_v131.py            # v0.13.1 welcome refs + eta guidance + soft bootstrap (~25 PASS)
   uv run tests/test_v132.py            # v0.13.2 welcome URL anti-pattern + validator note (~20 PASS, fast)
   uv run tests/test_v133.py            # v0.13.3 bootstrap_n default pinning (~5 PASS, fast)
   ```

   All six should end with "ALL TESTS PASSED" before any Claude Desktop setup.

## Running under Claude Desktop (manual install, dev mode)

The fastest iteration loop uses a direct mount in
`claude_desktop_config.json` rather than repacking a `.mcpb` after every
edit:

```json
{
  "mcpServers": {
    "brier": {
      "command": "uv",
      "args": ["run", "--directory", "/absolute/path/to/BRIER_MCP/mcp", "server.py"],
      "env": {
        "BRIER_RSCRIPT": "/absolute/path/to/Rscript"
      }
    }
  }
}
```

Quit Claude Desktop fully (tray icon -> Quit) and reopen between edits.

If you also have a `.mcpb`-installed version of the extension, the two
will conflict on the server name `brier`. Disable the extension while
developing via direct mount.

## Building a release bundle

```bash
mcpb validate manifest.json          # schema check
mcpb pack . brier-<version>.mcpb
```

The Anthropic `mcpb` CLI is npm-distributed:
`npm install -g @anthropic-ai/mcpb`. End users do not need it.

## Architecture invariants (don't change without re-testing)

These are commented in `server.py` but worth flagging here too:

1. `subprocess.run(..., stdin=subprocess.DEVNULL)` on every Rscript call.
   Without this, Rscript can inherit a parent stdin bound to the MCP
   stdio channel under Claude Desktop and stall on TTY probes.
2. Rscript flags are `--no-save --no-restore --no-init-file`, NOT
   `--vanilla`. `--vanilla` also implies `--no-environ` which suppresses
   `R_LIBS_USER` and breaks user-installed packages on Windows.
3. Rscript discovery falls back to well-known install locations because
   Claude Desktop subprocesses do not inherit the user's shell PATH on
   macOS.

## Naming convention

- **Technical identifiers** (manifest `name`, file names, env vars,
  Python module names, tool prefixes): lowercase `brier` /
  snake_case `brier_full`, `brier_i`, etc. Required by the MCPB schema.
- **User-visible strings** (manifest `display_name`, prose in
  docstrings, error messages, INSTALL.md): capitalized `BRIER`. This is
  what users see in Claude Desktop.

## Tool roadmap

| Version | Tools |
|---------|-------|
| v0.1.0 | `inspect_data` |
| v0.2.0 | + `list_data_directory`, `brier_i`, `brier_i_cv`, `brier_i_selection` |
| v0.3.0 | + `brier_predict`, `brier_evaluate` |
| v0.4.0 | + `brier_full`, `brier_full_selection` |
| v0.5.0 | + `brier_s`, `brier_s_selection`, `cal_ld`, `get_ldb` |
| v0.6.0 | + `start_analysis` (wizard, interview-style) |
| v0.7.0 | + `inspect_user_data`, data-first wizard reflow |
| v0.7.1 | M-rule, coarser BRIERfull eta grid, eta=0 baseline auto-fix, time_expectation, BRIERs un-standardize, output dir config |
| v0.8.0 | M=1 auto-ind, M-aware BRIERfull eta grid, nested-external detection, BRIERs boundary warn, y_train_expr stashing, deny-list safe namespaces |
| v0.8.1 | Revert y_train_expr / auto-stash / auto-apply (BRIERs un-stand is multi-option), wizard cross_family_comparison guidance, preprocessI/preprocessS wrappers |
| v0.8.2 | llms.txt corrections (upstream BRIER repo task) |
| v0.9 | Plot wrappers: brier_plot_eta with M=2 auto-heatmap, brier_plot_box, brier_plot_importance; PNG + CSV outputs |
| v0.10 | summarize_fit: comprehensive HTML report + standalone reproduce.R script |
| v0.10.1 | Wizard UX patch: structured multi-select options + render hints; directory-path handling |
| v0.10.2 | Multi-file `data_paths` across all data-loading tools; reproduce.R multi-file aware |
| v0.10.3 | Principled eta default + boundary diagnostic + `brier_plot_selection` (no-test-set eta plot) + per-call `output_dir` override |
| v0.11 | `brier_auto_tune_eta` for automatic eta-ceiling escalation across a fixed ladder |
| v0.12 | Auto-tune refinements: strict-equality escalation trigger + single-shot de-escalation |
| v0.13 | `prep_data` tool: 9 composable operations + audit log + summarize_fit integration |
| v0.13.1 | Welcome references restructure + eta_list anti-pattern guidance + soft bootstrap-plot degradation with fallback snippets |
| v0.13.2 | Wizard-text fixes: welcome URL-append anti-pattern called out by name; EXPRESSION VALIDATOR section explains the `::` whitelist (`::` was never blocked, but the assistant kept claiming so) |
| **v0.13.3** (current) | `bootstrap_n` default dropped from 100 to 20 in `summarize_fit`, `brier_plot_box`, `brier_plot_importance` for high-p exploratory workflows |
| v1.0.0 | Polish: INSTALL.md, GitHub release, signed `.mcpb` |
| v1.1.0 | SSH remote launch |
| v1.2.0 | PLINK, bcftools, GWAS Catalog, PGS Catalog wrappers |
| v1.3.0 | Skills layer for AI assistants |

## Cache directory

`brier_i` writes fit objects to `~/.cache/brier-mcp/fits/` (or
`$XDG_CACHE_HOME/brier-mcp/fits/` if set, or
`%LOCALAPPDATA%\brier-mcp\fits\` on Windows). Subsequent
`brier_i_selection` calls reload from there. Safe to delete at any
time; users will just need to refit.
