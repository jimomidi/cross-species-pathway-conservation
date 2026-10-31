# Cross-Species Pathway Conservation Toolkit

A modular Python toolkit for comparing somatic-mutation and functional
pathway-enrichment results between two species, and visualizing how
conserved the affected biological pathways are between them. It was built
for a two-species comparative-genomics workflow (a companion-animal
species compared with human data) but every script accepts its own inputs
via command-line arguments, and the species labels are centralized in
[`config.py`](config.py) so the toolkit can be adapted to a different
species pair.

## What this toolkit does

The pipeline covers four stages, each in its own folder. Every script is
independently runnable from the command line; none of them require a
notebook or IDE.

1. **Ortholog mapping** (`data_preparation/`) -- convert an Ensembl
   BioMart orthology export to CSV, map a query-species gene list to its
   orthologs in a reference species, and prepare deduplicated gene lists
   for functional-enrichment tools such as DAVID.
2. **Variant analysis** (`variant_analysis/`) -- compare the single-
   nucleotide substitution spectrum and the functional (exonic-consequence)
   mutation-type spectrum between species, using paired tumor/liquid-biopsy
   VCF files and ANNOVAR annotations for the query species.
3. **Pathway enrichment comparison** (`pathway_enrichment/`) -- compare
   DAVID functional-annotation results (GO Biological Process, KEGG,
   Reactome) between species, identify conserved and species-specific
   pathways, and test whether the observed conservation exceeds chance
   using a real pathway-universe size retrieved from public APIs.
4. **Figures** (`figures/`) -- publication-style plots built from the
   outputs of the previous stages (scatter, dumbbell, heatmap/lollipop,
   gene-overlap network maps, and a theme-grouped conceptual summary),
   plus a few standalone template figures.

A driver-gene cross-reference utility (`driver_gene_annotation.py`) and a
generic image-background utility (`utils/make_background_transparent.py`)
are also included.

Each stage's scripts read files and write results explicitly via
command-line arguments -- there is no hidden shared state between scripts.
Downstream figure scripts expect the specific output filenames produced by
`pathway_enrichment/compare_enrichment_results.py` (see that section
below).

## Repository layout

```
config.py                          Shared paths, species labels, and figure style.
driver_gene_annotation.py          Cross-reference genes against COSMIC CGC / IntOGen.
data_preparation/
    convert_biomart_export.py      Clean a raw BioMart export into CSV.
    map_orthologs.py               Map a gene list to orthologs and summarize the mapping.
    prepare_enrichment_gene_lists.py  Build gene lists for enrichment tools (3 subcommands).
variant_analysis/
    vcf_utils.py                   Shared VCF-parsing helpers.
    substitution_spectrum.py       Compare 12-class SNV substitution spectra.
    functional_mutation_types.py   Compare ANNOVAR functional mutation-type spectra.
pathway_enrichment/
    compare_enrichment_results.py  Compare DAVID results, find conserved pathways.
    pathway_universe_enrichment_tests.py  Fisher's exact / hypergeometric tests.
figures/
    plotting_utils.py              Shared plotting/IO helpers used by the figure scripts.
    cross_species_scatter_plot.py  Fold-enrichment scatter plot with correlation stats.
    cross_species_dumbbell_plot.py Per-pathway dumbbell comparison plot.
    pathway_network_map.py         Gene-overlap pathway similarity network (spring layout).
    grouped_circular_pathway_map.py  Same network, grouped circular/chord layout.
    heatmap_lollipop_plot.py       Fold-enrichment heatmap and log2-ratio lollipop plot.
    conceptual_theme_summary.py    Pathways grouped into biological themes, summary panel.
    ortholog_identity_scatter.py   Reciprocal ortholog sequence-identity scatter plot.
    evidence_legend.py             Standalone interaction-evidence legend (no data input).
    venn_diagram_template.py       Blank two-circle Venn diagram template (no data input).
utils/
    make_background_transparent.py  Batch-convert white image backgrounds to transparent.
data/                              Place your input files here (not tracked by git).
results/                           Generated tables and figures are written here (not tracked by git).
```

## Requirements

* Python 3.10+
* See [`requirements.txt`](requirements.txt) for package versions:
  pandas, numpy, matplotlib, matplotlib-venn, scipy, networkx, openpyxl,
  requests, Pillow.

Some figure scripts use "Times New Roman" as the default font
(`config.FONT_FAMILY`). If it is not installed on your system, matplotlib
will fall back to a default serif font; a few standalone figures use an
explicit serif font-family list with fallbacks already configured.

## Installation

```bash
git clone <this-repository-url>
cd cross-species-pathway-conservation
python -m venv .venv
source .venv/bin/activate    # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Input data

No example or sample data is included in this repository. Place your own
input files under `data/` (or point each script at any location via its
command-line arguments) in the formats described in each script's module
docstring:

| Stage | Script | Key inputs |
|---|---|---|
| Ortholog mapping | `convert_biomart_export.py` | Raw BioMart "Export all results to" text file |
| Ortholog mapping | `map_orthologs.py` | Query gene list (CSV) + BioMart CSV export |
| Ortholog mapping | `prepare_enrichment_gene_lists.py` | Mapped-ortholog table and/or a per-sample mutation table |
| Variant analysis | `substitution_spectrum.py` | Reference SNV catalog (CSV) + paired `<sample>_Tumor.vcf.gz` / `<sample>_cfDNA.vcf.gz` files |
| Variant analysis | `functional_mutation_types.py` | ANNOVAR `.exonic_variant_function` files + paired VCFs + reference functional-type counts CSV |
| Pathway enrichment | `compare_enrichment_results.py` | Per-species DAVID export CSVs (GO BP, KEGG, Reactome) |
| Pathway enrichment | `pathway_universe_enrichment_tests.py` | The combined summary CSV from the script above; optionally a local GOA GAF file |
| Driver genes | `driver_gene_annotation.py` | Candidate gene list + COSMIC CGC and IntOGen driver tables |
| Figures | `figures/*.py` | Outputs of the stages above (see each script's docstring) |

Run any script with `--help` to see its exact arguments and defaults.

## Usage

Each script can be run directly with `python <script>.py ...`. A typical
end-to-end sequence:

```bash
# 1. Ortholog mapping
python data_preparation/convert_biomart_export.py --input data/biomart/mart_export.txt
python data_preparation/map_orthologs.py \
    --gene-list data/query_gene_list.csv \
    --mapping data/biomart/mart_export.csv \
    --output-dir results/ortholog_mapping

# 2. Variant analysis
python variant_analysis/substitution_spectrum.py \
    --reference-mutations data/reference_mutations/reference_snvs.csv \
    --vcf-dir data/variants \
    --output-dir results/substitution_spectrum

python variant_analysis/functional_mutation_types.py \
    --annotation-dir data/annotation \
    --vcf-dir data/variants \
    --reference-functional-types data/reference_mutations/reference_functional_types.csv \
    --output-dir results/functional_mutation_types

# 3. Pathway enrichment comparison
python pathway_enrichment/compare_enrichment_results.py \
    --species-a-dir data/enrichment/canine \
    --species-b-dir data/enrichment/human \
    --output-dir results/pathway_conservation

python pathway_enrichment/pathway_universe_enrichment_tests.py reference-universe \
    --summary results/pathway_conservation/tables/Combined_Pathway_Conservation_Summary.csv \
    --goa-file data/reference_annotation/goa_human.gaf.gz \
    --output results/pathway_conservation/tables/reference_universe_enrichment_tests.csv

# 4. Figures (built from the pathway-conservation tables above)
python figures/cross_species_scatter_plot.py \
    --input-dir results/pathway_conservation/tables \
    --output-dir results/figures/cross_species_scatter_plot

python figures/cross_species_dumbbell_plot.py \
    --input-dir results/pathway_conservation/tables \
    --output-dir results/figures/cross_species_dumbbell_plot
```

See each script's `--help` output (or its module docstring) for the
complete set of inputs, outputs, and options -- including
`pathway_network_map.py`, `grouped_circular_pathway_map.py`,
`heatmap_lollipop_plot.py`, and `conceptual_theme_summary.py`, which all
read the same `results/pathway_conservation/tables` output directory.

## Outputs

Scripts write CSV/Excel tables and, for the figure scripts, image files
(PNG plus a subset of TIFF/PDF/SVG, depending on the script) to the
`--output-dir` you specify (default locations are under `results/`).
Nothing is written outside the directory you provide.

## Reproducibility

* All configurable parameters (species labels, KEGG organism codes,
  taxonomy IDs, color/marker schemes, minimum sample thresholds, sorting
  and layout options) are declared as named constants at the top of each
  script or centralized in `config.py` -- there are no unlabeled "magic
  numbers" controlling the analysis.
* All paths are relative to the repository root or explicitly provided via
  command-line arguments; no absolute or machine-specific paths are used.
* Sample identifiers are discovered automatically from input filenames
  (`<sample>_Tumor.vcf.gz` / `<sample>_cfDNA.vcf.gz`) rather than
  hardcoded, so the toolkit works unchanged with any sample set that
  follows this naming convention.
* No statistical procedure in this toolkit relies on random sampling. The
  one place randomness affects an output is the node layout of
  `figures/pathway_network_map.py` (a NetworkX spring layout), which uses
  a fixed seed (`LAYOUT_SEED` in that script) so the figure is identical
  across runs.

## Customization

* **Different species pair:** update `SPECIES_A`, `SPECIES_B`,
  `KEGG_ORGANISM_CODE`, and `NCBI_TAXONOMY_ID` in `config.py`.
* **Figure styling:** `DATABASE_COLORS`, `DATABASE_MARKERS`, and
  `FONT_FAMILY` in `config.py` are shared across the figure scripts.
* **Tissue/disease keyword flags:** `TISSUE_KEYWORDS` in
  `driver_gene_annotation.py` is a plain, editable list -- it is an
  annotation aid only and does not affect the driver classification logic.
* **Biological themes:** `THEME_KEYWORDS` in
  `figures/conceptual_theme_summary.py` defines the keyword groups used to
  bucket pathways into summary themes; adapt it to the biology of your own
  pathway set.
* **Pathway subset for the network map:** `pathway_network_map.py` accepts
  an optional `--pathway-allowlist` text file (one pathway name per line)
  to restrict the network to a chosen subset instead of plotting every
  conserved pathway.

## Known limitations / inputs you must supply yourself

* `variant_analysis/functional_mutation_types.py` does not compute the
  reference-species (species B) functional-mutation-type counts itself;
  this must be supplied as a two-column CSV via
  `--reference-functional-types`. Bring your own counts (for example,
  derived from a public annotated mutation catalog for the reference
  species).
* The same script's multi-transcript tie-breaking function
  (`choose_best`, used when an indel matches more than one candidate
  functional-consequence annotation) is a documented, explicit default
  (majority vote, with a fixed severity-order tiebreak) rather than a
  value carried over from an original analysis -- see the docstring in
  `functional_mutation_types.py` for the exact rule, and confirm it
  matches how you want multi-transcript annotations handled before relying
  on downstream counts.

## License

Released under the MIT License. See [`LICENSE`](LICENSE) for details.

Copyright (c) 2026 Javad Omidi
