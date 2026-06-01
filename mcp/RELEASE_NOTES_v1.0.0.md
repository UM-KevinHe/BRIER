# BRIER MCP - v1.0.0

First public-stable release of the BRIER MCP (Model Context Protocol) server. This extension lets Claude Desktop drive the [BRIER R package](https://github.com/UM-KevinHe/BRIER) for transfer-learning genetic risk prediction. Describe your data and your question in natural language; Claude inspects the file, picks the appropriate BRIER variant (`BRIERi` / `BRIERfull` / `BRIERs`), calls the corresponding R function on your machine, and explains the results.

**Your data never leaves your computer.** The server passes only file paths to R; raw data is read locally and only the resulting summaries (coefficients, validation metrics, plots) are returned to Claude.

---

> ### Before you start - three points that account for most setup issues
>
> 1. **You must install two prerequisites before this MCP works.** R (>= 4.0) and the BRIER R package. Python is NOT required; `uv` handles that. Details below.
>
> 2. **You configure this via `claude_desktop_config.json`, not via a one-click installer.** v1.0.0 ships as a git-clone-and-go server. A signed `.mcpb` bundle (the one-click drag-and-drop format) is planned for a later release once the feature surface stabilizes further.
>
> 3. **If Claude proposes the wrong analysis, force the wizard.** Say *"Please use the start_analysis wizard first."* That deterministically routes through the inspect-and-recommend path.

---

## What's in v1.0.0

- **27 MCP tools** covering: data inspection and wizard-guided routing; the three BRIER fitting families (`BRIERi`, `BRIERfull`, `BRIERs`) plus their selection/CV variants; automatic eta-grid escalation (`brier_auto_tune_eta`); prediction and evaluation; four plot types (eta curve, bootstrap comparison, variable importance, selection criterion); HTML reports with a runnable `reproduce.R` companion (`summarize_fit`); and a composable, leakage-aware data-prep pipeline (`prep_data`) with audit logging.
- **Data-first wizard.** `start_analysis` asks for a file path, inspects the file, derives outcome family / predictor type / sample sizes heuristically, and produces a tentative recommendation with the canonical call shape already filled in. The user confirms or corrects in plain English.
- **Honest cost warnings.** Bootstrap plots and high-`p` fits get heads-up notices about expected run time. If a plot exceeds a (user-configurable) timeout, the report still completes; the missing plot is replaced by a standalone R snippet the user can run themselves.
- **Cross-platform.** Tested install paths on macOS (Intel + Apple Silicon) and Linux; Windows install path documented but less heavily tested.
- **Local computation.** All fitting runs through the user's R installation. Only tool names, arguments, and result summaries transit the Claude Desktop API.

---

## Prerequisites

### 1. R (>= 4.0)

**Check whether R is already installed.**

- **macOS / Linux:** in a terminal, run `R --version`. A version number indicates R is installed.
- **Windows:** open the Start menu and type "R". If "R x64 4.x.x" or "RStudio" shows up, R is installed.

**If R is missing or older than 4.0**, install from <https://cran.r-project.org/>. Verify with `R.version.string`.

### 2. The `BRIER` R package

In R:

```r
install.packages("remotes")
remotes::install_github("UM-KevinHe/BRIER")
library(BRIER)   # should return to the prompt with no error
```

The first install takes a couple of minutes because BRIER compiles C++ code (via `RcppArmadillo`).

### 3. `uv`

The Python package/runtime manager used to launch the MCP server. Install from <https://docs.astral.sh/uv/getting-started/installation/>:

- **macOS:** `brew install uv` (or the official installer)
- **Linux:** the official installer
- **Windows:** `winget install astral-sh.uv` (or the official installer)

### What you do NOT need

- **Python.** `uv` manages a Python runtime for the server automatically.
- **Node.js** or **Conda / Anaconda.**
- **A C++ compiler.** `install.packages` offers Rtools (Windows) or Xcode Command Line Tools (macOS) on demand.

---

## Install the MCP server

```bash
git clone https://github.com/UM-KevinHe/BRIER.git
cd BRIER/mcp

# Optional: pin to this release
# git checkout v1.0.0

uv sync   # installs `mcp` and friends into an isolated venv
```

If you don't want git, download the source tarball from the Assets section below and `tar -xzf` it instead.

## Configure Claude Desktop

Edit your Claude Desktop config file:

- **macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows:** `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux:** `~/.config/Claude/claude_desktop_config.json`

Add a `brier` entry under `mcpServers`:

```json
{
  "mcpServers": {
    "brier": {
      "command": "uv",
      "args": [
        "run",
        "--directory",
        "/absolute/path/to/BRIER/mcp",
        "server.py"
      ],
      "env": {
        "BRIER_RSCRIPT": "/absolute/path/to/Rscript"
      }
    }
  }
}
```

**Fill in two values:**

- **`--directory`**: the absolute path of the `mcp/` folder you cloned (no `~`).
- **`BRIER_RSCRIPT`**: absolute path of your R installation's `Rscript`. Typical locations:
  - **macOS:** `/Library/Frameworks/R.framework/Resources/bin/Rscript` (or `/opt/homebrew/bin/Rscript` if installed via Homebrew)
  - **Linux:** `/usr/bin/Rscript` (or run `which Rscript`)
  - **Windows:** `C:\Program Files\R\R-4.x.x\bin\Rscript.exe` (substitute your version)

`BRIER_RSCRIPT` is optional; the server auto-discovers `Rscript` at standard locations if it's not set. Set it explicitly if you have multiple R versions installed.

## Restart Claude Desktop

**Fully quit and reopen.** New chats see the BRIER tools; existing chats won't.

---

## Verify

Open a new chat and say:

> Use BRIER to start an analysis on this file: `~/Desktop/test_brier.rds`

If you don't have a dataset to try, stage the synthetic example that ships with the BRIER R package:

```r
# In R
library(BRIER)
data(Data_BRIERi)
saveRDS(Data_BRIERi, "~/Desktop/test_brier.rds")
```

Claude should respond with the BRIER welcome message, ask one or two clarifying questions, then call `inspect_user_data` to peek at the file. If it does, the install is working.

---

## Troubleshooting

**Claude doesn't see BRIER tools after installing.**
Fully quit Claude Desktop (not just close the window) and reopen. Then open a *new* chat - existing chats won't pick up newly-added MCP servers.

**"Rscript not found" or "Could not find R".**
The `BRIER_RSCRIPT` env var is wrong, or it's unset and auto-discovery failed. Open `claude_desktop_config.json` and set `BRIER_RSCRIPT` to the absolute path of `Rscript` (use `which Rscript` on macOS/Linux to find it). On Windows the path must end in `Rscript.exe`.

**"package 'BRIER' is not available" or "There is no package called 'BRIER'".**
The BRIER R package isn't installed under the R version the MCP is calling. Open the R that `BRIER_RSCRIPT` points to, then run:

```r
remotes::install_github("UM-KevinHe/BRIER")
library(BRIER)
```

If you have multiple R versions installed, install BRIER under each one that you might configure the MCP against.

**Claude proposes the wrong tool.**
MCP tools are loaded lazily by Claude Desktop. Say *"Please use the start_analysis wizard first"* - that deterministically routes through the inspect-and-recommend path.

**A bootstrap plot times out on a high-`p` fit.**
`summarize_fit` in v1.0.0 defaults to `bootstrap_n=20` (down from 100 in earlier versions) precisely for this. If it still times out, ask Claude to retry with smaller `bootstrap_n` or to skip the bootstrap plots (`include_box_plot=False`, `include_importance_plot=False`). The eta and selection plots are always cheap.

**"R output was not valid JSON" errors out of nowhere.**
Check disk space. The fit cache (`~/.cache/brier-mcp/`) grows over time; if `/tmp` fills up, R can't create its working directory. Clear the cache with `rm -rf ~/.cache/brier-mcp/fits/*`.

**Any other issue.**
Open an issue at <https://github.com/UM-KevinHe/BRIER/issues> with your OS, R version (`R.version.string`), exact error text, and the prompt that triggered it.

---

## Privacy

- **Data files stay on your machine.** The server passes a file path to R, which reads the file locally.
- **Analysis results travel through Anthropic's API.** Coefficients, validation metrics, plot captions, and explanations appear in Claude's responses and are logged by Claude Desktop as part of the normal conversation.
- **Tool names and arguments transit Anthropic's API** as part of standard MCP operation.
- If your data file contains PHI or other sensitive information, plan accordingly: derived numerical results (e.g. *"coefficient for variable Age = 0.31"*) appear in Claude's responses.

---

## Uninstall

1. Remove the `brier` entry from `claude_desktop_config.json` and restart Claude Desktop.
2. Delete the cloned repo if desired: `rm -rf /path/to/BRIER`.
3. Optionally delete the cache: `rm -rf ~/.cache/brier-mcp/`.

The BRIER R package on your system remains installed and can be removed in R with `remove.packages("BRIER")`.

---

## Roadmap

- **v1.1.x**: SSH remote launch (drive a BRIER MCP running on a compute cluster from a local Claude Desktop). External-tool wrappers for PLINK, bcftools, GWAS Catalog, PGS Catalog. A skills layer for AI assistants.
- **v2.0.0**: Signed `.mcpb` bundle (one-click drag-and-drop install for Claude Desktop). Distribution polish.

Breaking changes follow semver: tool-signature changes, file-layout changes, or env-var changes will bump the major version. Additive features (new tools, optional kwargs) ship in minor versions.

---

## Full development history

Per-version dev notes (v0.1 through v0.13.3) are in [DEVLOG.md](DEVLOG.md). v1.0.0 is functionally identical to v0.13.3 plus a new top-level README, LICENSE, and version bump.

**Test suites passing in v1.0.0**: 21 (run `for t in tests/test_*.py; do uv run "$t"; done` from `mcp/`).

**Tool count**: 27.
