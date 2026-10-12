# DATA405 · Python & R for Genomic Data Science

Course website for Lecture 10 (Python and R), Lecture 11 (Probability and Statistics for Genomic Data) and Lab 1
of DATA405 Genomic Data Science, Howard University, Fall 2026 (Jinrui Xu).

All code cells run in the reader's browser via [quarto-live](https://github.com/r-wasm/quarto-live)
(Pyodide for Python, webR for R). No server, no installation.

## Layout

- `*.qmd` — pages (index, notes, rosetta, lab1, stats, data, setup)
- `data/` — the real, small datasets (see `data.qmd` for provenance)
- `downloads/` — slides (PDF), lab handout (PDF/DOCX), data zip
- `tools/build_data.py` — how the data package was built
- `docs/` — rendered site (GitHub Pages serves this folder)
- `_extensions/live/` — the quarto-live extension (vendored)

## Rebuild

Install [Quarto](https://quarto.org) ≥ 1.4, then `quarto render`. The output goes to `docs/`.

## Publish

GitHub → Settings → Pages → Source: *Deploy from a branch*, branch `main`, folder `/docs`.
