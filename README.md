# Inhibiting LINE-1

Manuscript and the calculations that draw its figures.

`paper/line1-inhibition/latex/inhibiting_line1.tex` is the two-column draft. `inhibiting_line1_article.tex` is the article layout. Both compile against the PNGs next to that directory. The PDFs in `latex/` are the current builds.

Figure scripts live in `analysis/nrti_line1_hiv/` and write into `paper/line1-inhibition/`:

- `plot_genomes.py` draws the LINE-1 and HIV genome figure.
- `plot_line1_hiv.py` fits LINE-1 versus HIV potency from `potency_data.tsv`.
- `plot_alimatravir_decay.py` draws intracellular triphosphate after an 11 mg dose.

Panel B of the nucleoside figure is the structure drawing `latex/alimatravir_stacked.svg`.

Regenerating the PNGs needs Python 3 with the packages in `requirements.txt`, and Google Chrome at `/Applications/Google Chrome.app`. From the repository root:

```bash
python analysis/nrti_line1_hiv/plot_genomes.py
python analysis/nrti_line1_hiv/plot_line1_hiv.py
python analysis/nrti_line1_hiv/plot_alimatravir_decay.py
```
