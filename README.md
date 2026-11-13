# Reproduce the research in the attached reference manuscript, Target Fishing and Molecular Docking of Terpenes. Recapitulate the 23-terpene, four-tool target-fishing consensus design and five-target docking panel (CB2/6KPF, androgen receptor/3V49, CYP19A1/3EQM, estrogen receptor beta/1L2J, and 11-bet

Terpedia study: `study_be6fa4a4-c9eb-45e6-916e-aec999992e4f`

Research snapshot SHA-256: `a0f2569040758f879017e98239620e7a280b83b1830c73f5d3a74545f5b7e325`
Notebook SHA-256: `798f59508252c76aadac39a8ea1b8aa1c81a86555057ed0da7d67035135f821e`

[Open the reproducibility notebook in Google Colab](https://colab.research.google.com/github/Terpedia/terpedia-reproduce-the-research-in-the-attached-reference-manusc-be6fa4a4/blob/main/notebooks/study_be6fa4a4-c9eb-45e6-916e-aec999992e4f-reproducibility-r11.ipynb)

## Reviewer informed revision 11

Revision 11 preserves the frozen revision 10 data and adds audits responding to both reviewer reports. It exports target-frequency and pair-level discrepancies, ambiguous target identifiers, tool-support counts without assuming missing tool outputs are negatives, target-specific docking labels, exploratory OLS statistics, terpene-only computational enrichment, and recorded receptor preparation. The notebook also offers structure drawings and PAINS screening when verified SMILES and RDKit are available.

Run all cells in Colab to regenerate the analysis. Frozen-data analyses require no API credentials. The optional natural-language API and private artifact retrieval require their own access. The notebook reconstructs recorded computations; it does not rerun docking or establish experimental binding.

The input includes 1,546 prediction rows and 260 author-consensus rows. Ambiguous identifiers, missing tool-attempt logs, original MOE settings, reference affinity provenance, ligand structures and redocking validation remain explicit revision tasks. See the final notebook section for the reviewer response checklist.

Generated CSV and JSON audit outputs are retained in `analysis/reviewer-audit/`. Rebuild the notebook with `python scripts/update_reviewer_notebook.py`. Revision 10 and earlier remain preserved.

Validation: all code cells executed locally against the frozen snapshot; notebook schema validated. The optional PAINS branch reported missing verified structures rather than generating results.
