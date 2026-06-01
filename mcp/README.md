# BRIER MCP

**An MCP (Model Context Protocol) server that lets Claude Desktop drive the [BRIER](https://github.com/UM-KevinHe/BRIER) R package for transfer-learning genetic risk prediction.** Describe your data and your question in natural language; Claude inspects the file, picks the right BRIER variant (`BRIERi`, `BRIERfull`, or `BRIERs`), calls the corresponding R function on your machine, and explains the results.

**Your data never leaves your computer.** The server passes only file paths to R; raw data is read locally and only summaries (coefficients, validation metrics, plots) are returned to Claude.

> **About BRIER.** BRIER is a regularized regression framework with an extra term that pulls coefficients toward external information. It improves prediction in a target cohort by borrowing structure from larger, related external sources (other-ancestry GWAS, pretrained model coefficients, or pooled individual-level data), while controlling for negative transfer via a tunable integration weight (`eta`). Background: Choi et al. 2020 ([Nat Protocols PRS tutorial](https://www.nature.com/articles/s41596-020-0353-1)). BRIER docs: <https://um-kevinhe.github.io/BRIER/>.

---

## Status

Version **v1.0.0** (first public-stable release). 27 tools, 21 test suites passing. See [DEVLOG.md](DEVLOG.md) for the full development history.

## What you need

Before installing, have these on your machine:

1. **R (>= 4.0)** from <https://cran.r-project.org/>. Verify in R: `R.version.string`.
2. **The BRIER R package.** In R:
   ```r
   install.packages("remotes")
   remotes::install_github("UM-KevinHe/BRIER")
   library(BRIER)   # should return to the prompt without error
   ```
3. **`uv`** (Python package/runtime manager) from <https://docs.astral.sh/uv/getting-started/installation/>. macOS users can `brew install uv`; Windows users can `winget install astral-sh.uv`; the official installer also works on all platforms.

What you do **not** need: Python (uv installs it), Node.js, Conda, or a C++ compiler. R's `install.packages` will offer to install Rtools (Windows) or Xcode Command Line Tools (macOS) if needed for compiling.

---

## Install

```bash
# Clone the BRIER repo and switch into the MCP subfolder
git clone https://github.com/UM-KevinHe/BRIER.git
cd BRIER/mcp

# Optional: pin to a specific release tag
# git checkout v1.0.0

# Install dependencies (just `mcp` itself; uv handles the rest)
uv sync
```

If you want a specific released version without cloning the whole repo, download the source tarball from the [Releases page](https://github.com/UM-KevinHe/BRIER/releases) and `tar -xzf` it.

## Point Claude Desktop at it

Edit your Claude Desktop config file:

- **macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows:** `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux:** `~/.config/Claude/claude_desktop_config.json`

Add an entry under `mcpServers` (create the key if it doesn't exist):

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

Two values to fill in:

- **`--directory`** points at the `BRIER/mcp` folder you cloned. Use the absolute path (no `~`).
- **`BRIER_RSCRIPT`** points at your R installation's `Rscript` executable.
  - macOS: `/Library/Frameworks/R.framework/Resources/bin/Rscript` (or `/opt/homebrew/bin/Rscript` if R came from Homebrew)
  - Linux: `/usr/bin/Rscript` (or `which Rscript` in a terminal)
  - Windows: `C:\Program Files\R\R-4.x.x\bin\Rscript.exe` (substitute your R version)

`BRIER_RSCRIPT` is optional; if you omit it, the server auto-discovers `Rscript` at standard install locations. Setting it explicitly is recommended if you have multiple R installations.

Fully quit and reopen Claude Desktop. The 27 BRIER tools will be available in any new chat.

---

## Verify

Open a new chat in Claude Desktop and say:

> Use BRIER to start a transfer-learning analysis on this file: `/path/to/your/data.rds`

Claude should respond with the welcome message, ask one or two clarifying questions, then call `inspect_user_data` to peek at the file. If it does, the install is working.

If you don't have a dataset to try yet, the BRIER R package ships with a small synthetic example you can stage in two commands:

```r
# In R
library(BRIER)
data(Data_BRIERi)
saveRDS(Data_BRIERi, "~/Desktop/test_brier.rds")
```

Then in Claude Desktop:

> Use BRIER to fit a model on `~/Desktop/test_brier.rds`

---

## Quickstart: a worked example

The most common path through the tools is: **inspect → start_analysis → fit → select → predict → summarize**. End to end, on the canonical `Data_BRIERi` example:

1. Say: "Use BRIER on `~/Desktop/test_brier.rds`"
2. Claude calls `inspect_user_data` to find out the file is a list with `target$train$X`, `target$train$y`, `beta.external` - recognized as a `BRIERi` setup.
3. Claude calls `start_analysis` (with the inspection results) and gets back a tentative recommendation: use `brier_i`, gaussian family, the canonical call shape is filled in for you.
4. Claude confirms with you, then runs `brier_i` to fit. The result is cached under a `fit_id`.
5. Claude runs `brier_i_selection` to pick `eta` and `lambda` using a validation set or IC.
6. Optionally: `brier_predict` on held-out data, `brier_plot_eta` to see the validation curve, and `summarize_fit` to get a single HTML report with everything.
7. For tricky cases where the optimum looks like it sits at the eta grid boundary, Claude can call `brier_auto_tune_eta` to widen the search automatically.

You don't have to know any of this in advance - Claude knows which tool to call when, and the wizard guides the early choices.

---

## What's in the box

27 MCP tools, grouped:

- **Wizard + inspection (4):** `start_analysis`, `inspect_data`, `inspect_user_data`, `list_data_directory`
- **Fit (5):** `brier_i`, `brier_i_cv`, `brier_full`, `brier_s`, `brier_auto_tune_eta` (auto-escalates the eta grid when the optimum hits the boundary)
- **Selection (3):** `brier_i_selection`, `brier_full_selection`, `brier_s_selection`
- **Prediction & evaluation (2):** `brier_predict`, `brier_evaluate`
- **Plots (4):** `brier_plot_eta`, `brier_plot_box`, `brier_plot_importance`, `brier_plot_selection`
- **Reporting (1):** `summarize_fit` (HTML report + standalone reproduce.R script)
- **Data prep (2):** `prep_data` (9 composable operations: rename columns, derive correlations, subset to common SNPs, harmonize alleles, etc.), `prep_data_log` (audit log)
- **LD utilities (2):** `cal_ld`, `get_ldb`
- **Legacy alignment (2):** `preprocess_i`, `preprocess_s`
- **I/O (2):** `set_output_directory`, `get_output_directory`

Detailed docstrings live in `server.py`. Test coverage is in `tests/test_v*.py`.

---

## Configuration

### Where outputs go

By default, plots and reports land in a per-session cache directory (`~/.cache/brier-mcp/`). To control where they land, either:

- Tell Claude something like "save outputs to `~/Desktop/MyProject/`" - Claude will call `set_output_directory`.
- Pass `output_dir="/absolute/path"` to individual tools (`summarize_fit`, `brier_predict`, etc.).

### Cache management

The fit cache (`~/.cache/brier-mcp/fits/`) grows as you run more analyses. Several GB after a few dozen fits is normal. Clear it with:

```bash
rm -rf ~/.cache/brier-mcp/fits/*
```

If you ever see "R output was not valid JSON" errors out of nowhere, check disk space; the most common cause is a full `/tmp` because the cache grew large.

---

## Troubleshooting

**"BRIER package is not installed" or "package 'BRIER' is not available".**
The R package is missing or installed under a different R version than the one `BRIER_RSCRIPT` points at. Open the R version that `BRIER_RSCRIPT` points to and run `remotes::install_github("UM-KevinHe/BRIER")`.

**"Rscript not found".**
Either set `BRIER_RSCRIPT` to the full path of your Rscript executable in `claude_desktop_config.json`, or check that `Rscript` is on your PATH (`which Rscript` in a terminal).

**Claude doesn't see BRIER tools.**
Fully quit Claude Desktop (not just close the window) and reopen. New chats see the tools; existing chats may not.

**The assistant proposes the wrong analysis.**
Say *"Use the start_analysis wizard first"*. That deterministically routes through the inspect-and-recommend path rather than guessing.

**The fit times out on large data (`p >= 10,000`).**
Bootstrap plots in `summarize_fit` default to `bootstrap_n=20` in v1.0.0+ (down from 100 in earlier versions) for this exact reason. If a plot still times out, you can: (1) pass a smaller `bootstrap_n` explicitly, or (2) omit the box/importance plots (`include_box_plot=False`, `include_importance_plot=False`). The eta and selection plots are cheap regardless of `p`.

**For other issues**, open an issue at <https://github.com/UM-KevinHe/BRIER/issues> with: your OS, your R version (`R.version.string`), the exact error text, and the prompt that triggered it.

---

## Privacy

- **Data files stay on your machine.** The server passes file paths to R; R reads the file locally.
- **Analysis results travel through Anthropic's API** as part of the normal conversation. Coefficients, validation metrics, and plot captions appear in Claude's responses and are therefore logged by Claude Desktop.
- **Tool names and arguments transit Anthropic's API** as part of standard MCP operation.
- If your data file contains PHI or other sensitive information, plan accordingly: derived results (e.g. *"coefficient for variable Age = 0.31"*) will appear in Claude's responses.

---

## Development

To run the test suite:

```bash
cd BRIER/mcp
uv run tests/test_v133.py   # most recent specific suite; ~5 PASS, fast
# or run all:
for t in tests/test_*.py; do uv run "$t"; done
```

See [DEVLOG.md](DEVLOG.md) for the per-version development history (v0.1 through v1.0.0).

## License

This MCP wrapper is MIT-licensed (see [LICENSE](LICENSE)). The underlying BRIER R package has its own license; check its repository for terms.

## Citation

If you use BRIER in published work, please cite the BRIER paper (see <https://github.com/UM-KevinHe/BRIER>).
