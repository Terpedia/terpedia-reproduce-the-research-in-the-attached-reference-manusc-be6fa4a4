"""Build revision 11 from the preserved revision 10 snapshot."""
import copy, hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = next((root / 'notebooks').glob('*-r10.ipynb'))
nb = json.loads(source.read_text())
def md(text):
    return {'cell_type': 'markdown', 'metadata': {}, 'source': text.splitlines(True)}
def code(text):
    return {'cell_type': 'code', 'metadata': {}, 'source': text.splitlines(True), 'execution_count': None, 'outputs': []}
nb['cells'][0] = md('''# Reviewer informed audit of terpene target fishing and docking

Revision 11 extends revision 10 in response to the two reviewer reports supplied on 13 November 2026. The original frozen study snapshot and its checksum are preserved. Run all cells to reconstruct recorded results and generate the new audit tables. No new docking or binding experiments are performed by this notebook.

The panel contains 23 terpenes and five selected targets. The results describe computational method agreement within this panel; they do not establish binding affinity, biological activity, target selectivity, or general performance across terpenes. The 130 recorded Vina observations use a single seed, exhaustiveness 4 and three requested poses. Completed runs do not establish successful redocking validation.
''')
# Avoid implying experimental accuracy when comparing two computational label sets.
s = ''.join(nb['cells'][20]['source']).replace('precision_vs_author', 'fraction_recomputed_matching_author').replace('recall_vs_author', 'fraction_author_recovered').replace('author_recall', 'fraction_author_recovered')
s = s.replace('status="reproduced"', 'status="adapted"')
nb['cells'][20]['source'] = s.splitlines(True)
nb['cells'][24] = md('''### Historical enrichment reconstruction

The following curves reproduce the original convention that pools known ligands and TF-positive terpenes as positives. They measure agreement between computational labels and docking rank, not predictive accuracy or experimental activity enrichment. Pooling controls with predicted positives can alter the result; the reviewer audit below adds a terpene-only sensitivity analysis. Enrichment cutoffs use the historical row ordering to resolve tied scores; inspect sensitivity to ties before interpreting a cutoff. A missing prediction row is not proof that a tool was run and returned a negative result.
''')
nb['cells'][22]['source'] = ''.join(nb['cells'][22]['source']).replace('\"cyp19a1\": \"CYP19A1\",', '\"cyp19a1\": \"CYP19A1\", \"cyp19a1aromatase\": \"CYP19A1\",').splitlines(True)
nb['cells'][21] = md('''### Target-specific labels and molecular-size effects

The following historical reconstruction uses membership in the author consensus for each ligand–target pair. Absence from that list is a computational label, not experimentally established inactivity. The quantity called ligand efficiency below is absolute docking score per heavy atom; it is a descriptive score normalization, not measured thermodynamic ligand efficiency. Dividing by heavy-atom count does not demonstrate removal of size bias. MOE and Vina scores must remain on separate scales. A negative sign or nominal energy unit does not make a docking score an experimental binding free energy.
''')
new = []
new.append(md('''## Reviewer audit and revision priorities

Reviewer 1 raises ten concerns; Reviewer 2 expands the requirements for selection, consensus rules, receptor structures, refinement, statistics, tables, PAINS and CYP19A1. This section audits available records without changing the frozen observations. Items described as unresolved require source records or additional computation before the manuscript can claim they are resolved.

The tool-count rule here requires at least two distinct supporting tools. Consequently 2/2 and 2/4 both pass that rule, but they do not document identical tool coverage. The raw supplement lists positive predictions: it cannot establish whether an absent tool was unrun, failed, returned no targets or fell below a threshold. Tool outputs may also share reference data, so agreement is not evidence of statistical independence.
'''))
new.append(code('''from importlib.metadata import version, PackageNotFoundError
from datetime import datetime, timezone
import numpy as np
from scipy import stats

review_output = analysis_output / "reviewer-audit"
review_output.mkdir(parents=True, exist_ok=True)
def audit_table(frame, filename):
    frame.to_csv(review_output / filename, index=False)
    display(frame)
    return frame

# Preserve all source rows. Keys only reconcile identifiers already in this snapshot.
identity_rows = pd.concat([target_predictions, author_consensus], ignore_index=True)
alias_candidates = {}
for row in identity_rows.to_dict("records"):
    identifiers = [str(row[k]) for k in ("chembl_id", "uniprot_id") if pd.notna(row.get(k)) and str(row[k]).strip()]
    if identifiers:
        alias_candidates.setdefault(target_key(row.get("target_name")), set()).add(identifiers[0])
ambiguous_aliases = {k: sorted(v) for k, v in alias_candidates.items() if len(v) > 1}
def audit_target_id(row):
    for field in ("chembl_id", "uniprot_id"):
        if pd.notna(row.get(field)) and str(row[field]).strip():
            return str(row[field])
    name = target_key(row.get("target_name"))
    candidates = alias_candidates.get(name, set())
    return next(iter(candidates)) if len(candidates) == 1 else "name:" + name

raw = target_predictions.copy()
auth = author_consensus.copy()
for frame in (raw, auth):
    frame["audit_target_id"] = frame.apply(audit_target_id, axis=1)
    frame["ligand_key"] = frame["terpene"].map(normalized_name)
raw_pairs = raw.groupby(["ligand_key", "audit_target_id"])["tool"].agg(lambda x: sorted(set(x))).reset_index(name="supporting_tools")
raw_pairs["support_count"] = raw_pairs["supporting_tools"].map(len)
raw_pairs["passes_two_tool_rule"] = raw_pairs["support_count"] >= 2
raw_pairs["observed_tools_for_ligand"] = raw_pairs["ligand_key"].map(raw.groupby("ligand_key")["tool"].nunique())
raw_pairs["attempted_tools_denominator"] = pd.NA
raw_pairs["coverage_note"] = "Observed positive-output tools only; attempted-tool denominator unknown"
audit_table(raw_pairs, "pair-tool-support.csv")

A = set(map(tuple, auth[["ligand_key", "audit_target_id"]].drop_duplicates().to_numpy()))
R = set(map(tuple, raw_pairs.loc[raw_pairs.passes_two_tool_rule, ["ligand_key", "audit_target_id"]].to_numpy()))
pair_discrepancies = pd.DataFrame([
    {"ligand_key": ligand, "target_id": target, "comparison": "author_only" if (ligand, target) in A else "recomputed_only"}
    for ligand, target in sorted(A ^ R)
], columns=["ligand_key", "target_id", "comparison"])
audit_table(pair_discrepancies, "consensus-pair-discrepancies.csv")
frequency_check = target_frequencies.copy()
frequency_check["audit_target_id"] = frequency_check.apply(audit_target_id, axis=1)
author_counts = auth.drop_duplicates(["ligand_key", "audit_target_id"]).groupby("audit_target_id").size()
recomputed_counts = raw_pairs[raw_pairs.passes_two_tool_rule].groupby("audit_target_id").size()
frequency_check["author_unique_ligand_count"] = frequency_check.audit_target_id.map(author_counts).fillna(0).astype(int)
frequency_check["recomputed_two_tool_count"] = frequency_check.audit_target_id.map(recomputed_counts).fillna(0).astype(int)
unresolved_identity = frequency_check.audit_target_id.str.startswith("name:")
frequency_check.loc[unresolved_identity, ["author_unique_ligand_count", "recomputed_two_tool_count"]] = np.nan
frequency_check["identity_status"] = np.where(unresolved_identity, "unresolved_identity_do_not_compare_counts", "identifier_matched_pending_source_curation")
frequency_check["reported_minus_author"] = frequency_check.terpenes_predicting - frequency_check.author_unique_ligand_count
frequency_check["reported_minus_recomputed"] = frequency_check.terpenes_predicting - frequency_check.recomputed_two_tool_count
audit_table(frequency_check, "target-frequency-reconciliation.csv")
audit_table(pd.DataFrame([{"normalized_target_name": k, "candidate_ids": "; ".join(v)} for k, v in sorted(ambiguous_aliases.items())]), "ambiguous-target-identifiers.csv")
print("Ambiguous name-to-identifier mappings requiring curation:", len(ambiguous_aliases))
print("Recomputed zero counts mean no qualifying pair in supplied rows; missing outputs remain unresolved.")
'''))
new.append(md('''### Pair labels and statistical reporting

The audit table preserves author-consensus membership and raw tool support separately. Known ligands are controls, not TF-positive terpenes. No missing prediction is silently promoted to an experimentally inactive compound.

Ordinary least squares below fits absolute docking score as a function of heavy-atom count, separately for each method and target, with terpenes alone and controls pooled, plus the two terpene classes separately. It reports sample size, compounds, slope, signed Pearson r, R², slope confidence interval and an unadjusted two-sided slope p-value. These exploratory estimates assume independent observations and linear-model errors; small, selected panels and shared scaffolds limit inference. No confirmatory significance claims or cross-method score comparisons are made.
'''))
new.append(code('''support_lookup = {(r.ligand_key, r.audit_target_id): r.support_count for r in raw_pairs.itertuples()}
target_id_candidates = {}
for r in auth.itertuples():
    target_id_candidates.setdefault(target_key(r.target_name), set()).add(r.audit_target_id)
target_id_lookup = {k: next(iter(v)) for k, v in target_id_candidates.items() if len(v) == 1}
review_docking = docking.copy()
review_docking["author_tf_status"] = [
    "control_not_applicable" if r.ligand_class == "known_ligand" else
    ("TF+_author_consensus" if r.tf_positive else "absent_from_author_consensus")
    for r in review_docking.itertuples()
]
review_docking["raw_support_tools"] = [
    support_lookup.get((normalized_name(r.ligand), target_id_lookup.get(target_key(r.target))), pd.NA)
    if r.ligand_class != "known_ligand" else pd.NA for r in review_docking.itertuples()
]
review_docking["score_per_heavy_atom"] = review_docking.ligand_efficiency
review_docking["experimental_affinity_value"] = pd.NA
review_docking["experimental_affinity_source"] = "Requires curated affinity type, unit, assay, species and citation"
audit_table(review_docking, "target-specific-docking-table.csv")

regression_rows = []
for (method, target), group in docking.groupby(["method", "target"]):
    for label, subset in [
        ("all_with_controls", group),
        ("terpenes_only", group[group.ligand_class != "known_ligand"]),
        ("monoterpenes_only", group[group.ligand_class == "monoterpene"]),
        ("sesquiterpenes_only", group[group.ligand_class == "sesquiterpene"]),
    ]:
        subset = subset.dropna(subset=["heavy_atoms", "score_magnitude"])
        row = {"method": method, "target": target, "group": label, "n": len(subset),
               "compounds": "; ".join(subset.ligand), "x": "heavy_atoms", "y": "absolute_docking_score"}
        if len(subset) >= 3 and subset.heavy_atoms.nunique() > 1 and subset.score_magnitude.nunique() > 1:
            fit = stats.linregress(subset.heavy_atoms, subset.score_magnitude)
            margin = stats.t.ppf(.975, len(subset)-2) * fit.stderr
            row.update(slope=fit.slope, intercept=fit.intercept, pearson_r=fit.rvalue,
                       r_squared=fit.rvalue**2, slope_ci95_low=fit.slope-margin,
                       slope_ci95_high=fit.slope+margin, p_unadjusted=fit.pvalue,
                       status="exploratory_ols")
        else:
            row["status"] = "undefined_insufficient_n_or_constant_variable"
        regression_rows.append(row)
regressions = audit_table(pd.DataFrame(regression_rows), "size-score-regressions.csv")

# Sensitivity: remove known controls; computational labels remain author membership only.
sensitivity = []
for (method, target), group in docking.groupby(["method", "target"]):
    population = group[group.ligand_class != "known_ligand"]
    for row in enrichment_curve(population):
        sensitivity.append({"method": method, "target": target, "population": "terpenes_only",
                            "label_source": "author_consensus_membership", **row})
audit_table(pd.DataFrame(sensitivity), "terpene-only-computational-enrichment.csv")
'''))
new.append(md('''### Receptor and validation evidence

The source PDB structures and any recorded preparation parameters are exported below. Co-crystal-derived boxes depend on a ligand-stabilized receptor conformation; this workflow is not agnostic to known ligands. Neither run completion nor the presence of a co-crystal control establishes successful pose recovery. Inspect validated redocking outputs, symmetry-aware heavy-atom RMSD, pose rank, atom mapping and coordinates before claiming validation.

For each target, the manuscript still needs structure-selection justification, species, resolution, mutations, site completeness, waters and cofactor treatment. Original MOE Triangle Matcher, London dG, GBVI/WSA dG, induced-fit and any MM/PBSA steps require their actual chronological settings and final-score definition. They cannot be inferred from the Vina workflow. CYP19A1 additionally requires documented heme/iron parameterization and coordination treatment; docking here does not establish selectivity against other CYP isoforms.
'''))
new.append(code('''receptor_records = []
for run in STUDY["runs"]:
    outputs = run.get("outputs", {})
    if outputs.get("receptor"):
        receptor_records.append({"run_id": run.get("id"), "run_status": run.get("status"),
                                 "receptor": outputs["receptor"], "method": outputs.get("method", {}),
                                 "validation_assessment": "Requires explicit pose/RMSD evidence; completion alone is insufficient"})
receptor_audit = pd.json_normalize(receptor_records, sep=".")
audit_table(receptor_audit, "recorded-receptor-preparation.csv")
(review_output / "recorded-receptor-preparation.json").write_text(json.dumps(receptor_records, indent=2))
'''))
new.append(md('''### Structures and PAINS screening input

Provide `ligand-structures.csv` in the working directory with `ligand`, `smiles`, and `identity_source` columns after checking stereochemistry against the source records. The optional cell below draws available structures and evaluates RDKit PAINS substructure alerts. It records missing identities rather than resolving names automatically. It does not download or install software. In Colab, install RDKit first if needed (`%pip install rdkit`), then rerun the cell.

PAINS alerts cannot distinguish genuine biological promiscuity from interference on their own. Neither an alert nor its absence establishes activity or assay cleanliness. Report alerts as screening flags and plan relevant experimental counterscreens. Reference: [Baell and Nissink, PAINS utility and limitations](https://doi.org/10.1021/acschembio.7b00903).
'''))
new.append(code('''structures_path = Path("ligand-structures.csv")
structure_rows = pd.DataFrame({"ligand": sorted(docking.ligand.unique())})
if structures_path.exists():
    supplied = pd.read_csv(structures_path).fillna("")
    assert {"ligand", "smiles", "identity_source"} <= set(supplied), "Missing structure provenance columns"
    assert not supplied.ligand.duplicated().any(), "Resolve duplicate ligand identities before screening"
    structure_rows = structure_rows.merge(supplied[["ligand", "smiles", "identity_source"]], how="left", on="ligand").fillna("")
else:
    structure_rows["smiles"] = ""
    structure_rows["identity_source"] = ""
structure_rows["pains_status"] = "blocked_missing_verified_structure"
structure_rows["pains_alerts"] = ""
try:
    from rdkit import Chem, rdBase
    from rdkit.Chem import FilterCatalog, Draw
    parameters = FilterCatalog.FilterCatalogParams()
    parameters.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(parameters)
    molecules, names = [], []
    for idx, row in structure_rows.iterrows():
        if not row.smiles or not row.identity_source:
            continue
        molecule = Chem.MolFromSmiles(row.smiles)
        if molecule is None:
            structure_rows.loc[idx, "pains_status"] = "invalid_smiles"
            continue
        alerts = [match.GetDescription() for match in catalog.GetMatches(molecule)]
        structure_rows.loc[idx, "pains_status"] = "alert" if alerts else "no_alert_not_experimental_clearance"
        structure_rows.loc[idx, "pains_alerts"] = "; ".join(alerts)
        molecules.append(molecule)
        names.append(row.ligand)
    if molecules:
        drawing = Draw.MolsToGridImage(molecules, legends=names, molsPerRow=4, subImgSize=(250, 200))
        display(drawing)
        drawing.save(str(review_output / "ligand-structures.png"))
    print("RDKit version:", rdBase.rdkitVersion)
except ImportError:
    structure_rows.loc[structure_rows.smiles.ne(""), "pains_status"] = "blocked_rdkit_unavailable"
audit_table(structure_rows, "ligand-structure-pains-audit.csv")
'''))
new.append(md('''### Work remaining for the manuscript revision

| Reviewer request | Notebook contribution | Remaining requirement |
|---|---|---|
| R1 1 and 9; R2 consensus | Pair discrepancies, tool support and target-frequency reconciliation | Reconcile manuscript values against these source rows; curate aliases and actual tool run/no-result logs and thresholds |
| R1 2; R2 docking | Export recorded receptor and engine settings | Supply or perform validated redocking with poses, RMSD and atom mapping; document failed cases |
| R1 3; R2 methods | Preserve snapshot, parameters and preparation provenance | Recover complete original MOE/refinement settings and software citations; receptor selection rationale |
| R1 4 and 6 | Separate score scales; describe per-atom normalization cautiously | Rewrite affinity and efficiency claims; inspect remaining size dependence |
| R1 5; R2 controls | Explicit target-specific controls in docking table | Curate experimental affinities, assay context, species and control-selection criteria |
| R1 7 | Computational enrichment with terpene-only sensitivity | Avoid labeling computational agreement as predictive accuracy |
| R1 8; R2 selection | Report n and selected-panel limits | Document original collection, inclusion/exclusion, stereochemistry and selection rationale; justify five targets |
| R1 10 | Hypothesis-generation scope stated | Revise manuscript conclusions and specify biochemical/binding follow-up |
| R2 statistics | OLS slopes, signed r, R², CI, n and exploratory p-values | Match exact manuscript regressions; justify models and any inferential claims |
| R2 tables and PAINS | Target-specific labels and optional structure/PAINS workflow | Supply verified structures and review generated tables; explain screening limitations |
| R2 CYP19A1 | Export recorded cofactor treatment and explicit limits | Verify iron/coordination model and discuss CYP selectivity without asserting it |

The embedded supplementary data supersede the older local reproduction's statement that the full prediction supplement was unavailable. Presence of these rows does not establish completeness of tool attempts or resolve disagreements automatically. Save the generated audit folder with the revised manuscript and use its discrepancies to guide a point-by-point response.
'''))
new.append(code('''packages = {}
for package in ("numpy", "pandas", "scipy", "matplotlib", "rdkit"):
    try:
        packages[package] = version(package)
    except PackageNotFoundError:
        packages[package] = "not_installed"
review_manifest = {
    "notebook_revision": 11, "source_revision": 10,
    "review_basis": "Two reviewer reports supplied 2026-11-13",
    "source_snapshot_sha256": EXPECTED_RESEARCH_SHA256,
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "packages": packages,
    "new_docking_performed": False,
    "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(review_output.iterdir()) if p.is_file() and p.name != "audit-manifest.json"},
}
(review_output / "audit-manifest.json").write_text(json.dumps(review_manifest, indent=2))
print("Reviewer audit exported to", review_output)
'''))
nb['cells'].extend(new)
for index, cell in enumerate(nb['cells']):
    cell['id'] = 'r11-' + str(index)
    if cell['cell_type'] == 'code':
        cell['execution_count'] = None
        cell['outputs'] = []
nb['metadata']['colab']['name'] = 'Terpene target fishing and docking reviewer audit r11.ipynb'
nb['metadata']['terpedia']['reviewRevision'] = 11
nb['metadata']['terpedia']['sourceNotebookRevision'] = 10
out = source.with_name(source.name.replace('-r10.ipynb', '-r11.ipynb'))
out.write_text(json.dumps(nb, indent=2, ensure_ascii=False) + '\n')
print(out)
