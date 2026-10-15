# Supplementary Material

**Semiempirical exploration of Eu(III) coordination in microsolvated homogalacturonan fragments: effects of the methylesterification pattern**

Lucas Faria da Silva. Repository: https://github.com/lucasfffaria/euIII-homogalacturonan-coordination

This folder is the audit and reproducibility package of the article (revised version). It holds what is needed to check the results and nothing else: the MOPAC inputs and outputs of every calculation, the construction scripts and SMILES of the molecular models, the processed data behind every table and figure, and the scripts that connect one to the other. With it an independent reader can

1. inspect the input and the output of every optimization and single point of the article (`mopac_optimizations/`, `revision/`, `mopac_single_points/`);
2. confirm that the archived geometries and energies are those of the MOPAC outputs, and rebuild the descriptor table from them (`code/reproduce_descriptors_v2.py`);
3. recompute Tables 2–6, the supplementary tables and every number quoted in the text from the descriptor table (`code/reproduce_tables_v2.py`);
4. repeat the checks of the Eu–O pair term, of the vibrational analysis, of the displacement test and of the crystal-structure benchmark (Section 4);
5. regenerate the molecular models and their SMILES (`code/galA.py`, `assemble.py`, `make_campaign.py`, `smiles_table.py`).

Every file has one location. The only repeated contents are the per-job parameter files that MOPAC reads from the job folder (`eu_fix.txt`). What is deliberately not in the package is listed in Section 7.

---

## 1. Quick check (Python only, about one minute)

```bash
cd code
python reproduce_descriptors_v2.py --compare      # MOPAC outputs -> archived geometries -> descriptor table and Table S9
python reproduce_tables_v2.py --out tables_check  # descriptor table -> Tables 2-6, S2-S5, key numbers
```

The first command must end with `identical values = True`. The second writes one CSV per table; they must equal the files of `data/tables_v2/`.

---

## 2. Models and names

| Key in file and job names | Model (name in the article) | First keywords of the MOPAC input |
|---|---|---|
| `PM7` | Sparkle/PM7+EuO: Sparkle/PM7 with an explicit Eu–O pair term | `PM7 SPARKLE … EXTERNAL=eu_fix.txt` (`PM7 … EXTERNAL=eu_fix.txt` in the original campaign) |
| `SRM1` | Sparkle/RM1, as distributed | `RM1 SPARKLE` |
| `RM1` | RM1 lanthanide model (Eu with a 5d/6s/6p basis) | `RM1` |
| `PM7std` | Sparkle/PM7, as distributed (comparator) | `PM7 SPARKLE` |

`eu_fix.txt` contains the two lines `XFAC_Eu O 2.207500` and `ALPB_Eu O 2.500000`. Under PM7 the keyword `SPARKLE` does not change the calculation for Eu (identical heats of formation, Table S8); under RM1 it selects a different model. The remaining keywords of an optimization are `CHARGE=q AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET`.

Job names are `<model>_<pattern>_<conformer>_<start>[_cw0|_cw4]` (e.g. `SRM1_B3_h2_bid`); controls are `<model>_CTRL_<system>_<variant>` (e.g. `PM7_CTRL_Eu_w9` = [Eu(H₂O)₉]³⁺); structures of the crystal benchmark are `<model>_XBnn_<label>`.

| Code | Meaning |
|---|---|
| pattern | A6, H3blk, H3alt, P3, B3, L3, B2, D2, S1, M6 (Figure 1); column `pattern` gives the sequence (A = –COO⁻, M = –COOCH₃, H = –COOH, non-reducing → reducing end) |
| conformer | `h3` 3₁-like helix, `h2` 2₁-like helix, `cl` irregular coil |
| start (column `mode`) | `bid` bidentate carboxylate; `mono2` monodentate plus bridge; `hyd` hydrated ion away from the chain; `ester2` ester carbonyl of the fully esterified chain (M6); `ester` first version of that start, superseded (kept and flagged, not analysed); controls: `w8`–`w12`, `bid`, `hyd`, `ester` |
| `_cw0`, `_cw4` | hydration-sensitivity series: 0 or 4 hydrogen-bonded water molecules on each free carboxylate not bound to Eu (main series: 2) |
| `series` | `main`, `sens_hydration`, `control` |

---

## 3. Contents

### Documents

| File | Purpose |
|---|---|
| `README.md` | this file |
| `Supplementary_Information.docx` | Supplementary Information of the article: Sections S1–S5, Tables S1–S15, Figures S1–S9 |

### `mopac_optimizations/` – original campaign (Sections 2.4–2.5)

180 optimizations: 101 with Sparkle/PM7+EuO and 79 with the RM1 lanthanide model (main series, hydration series and controls). One folder per job.

| File | Purpose |
|---|---|
| `<job>.mop` | MOPAC input: keywords and starting geometry in mixed Cartesian/internal coordinates (O–H lengths as fixed internal coordinates) |
| `<job>_start.xyz` | starting geometry in MOPAC atom order; reference of the integrity criterion and of `chain_rmsd`, `max_dtorsion` |
| `<job>.out` | MOPAC output (status, final geometry, heat of formation, gradient, charges) |
| `<job>.arc` | MOPAC archive of the final structure (written for runs that ended normally) |
| `eu_fix.txt` | Eu–O pair parameters read through `EXTERNAL=` (Sparkle/PM7+EuO folders only) |
| `<job>_rerun.mop`, `.out`, `.arc` | repetition of a run that did not reach self-consistence, from the same starting geometry with `SHIFT=15 PULAY CAMP ITRY=999` added (23 jobs) |
| `manifest.json` | construction record of the 180 jobs (model, charge, start, target residues and distances, water molecules, residue atom map); read by `smiles_table.py` |
| `progress_campaign6.log` | log of the production run |

### `revision/` – calculations added in the revision

One sub-folder per stage; each holds a manifest (`manifest_R*.csv`: job, keywords, origin), the log of the runner (`progress.log`) and one folder per job with the same files as above.

| Stage | Content | Article |
|---|---|---|
| `R1_param_path` | 17 single points that identify which Eu–O pair parameters MOPAC uses and the effect of the keyword `SPARKLE` | Sections 2.3, S1, S3; Table S8 |
| `R2_sparkle_pm7_std` | Sparkle/PM7 as distributed: 98 optimizations from the starting structures of the original campaign | Section 3.2, Figure 7, Tables 2, S1, S2, S4 |
| `R4_sparkle_rm1_scans`, `R4b_sparkle_rm1_scans` | rigid Eu–O scans with Sparkle/RM1: Eu³⁺ and the six rigid ligands (7 single points) and the 132 scan points | Figure 3, Tables 2, S3c |
| `R5_sparkle_rm1` | Sparkle/RM1: 98 optimizations from the same starting structures | Sections 3.3–3.10, Tables 3–6, S1 |
| `R6_M6_ester_restart` | corrected ester-carbonyl starts of the fully esterified chain (`ester2`), four models, 12 optimizations, with their `_start.xyz` | Sections 2.4, 3.8; Figure 5 |
| `R3_force`, `R8_force_new_series` | numerical Hessians (`FORCE`) at the final geometries of the structures analysed: inputs of all 230 calculations; outputs (`.out`, `.aux`) of nine representative ones (Section 7, item 2) | Sections 2.5, 3.11, S4; Figure 17a; Table S10 |
| `R7_displace_reopt`, `R9_displace_reopt_2` | displacement test: 186 re-optimizations after a displacement of 0.25 Å along the imaginary modes | Section 3.11, S4; Figure 17b,c; Table S11 |
| `R10_xray_benchmark` | crystal-structure benchmark: 17 complexes × 4 models (68 optimizations); `_reference/` holds the crystal coordinates of the isolated complexes (`*_xray.xyz`), the table of structures (`R10_structures.csv`: COD identifier, DOI, charge, ligands as SMILES, donors) and the 185 reference Eu–O distances (`xray_EuO_distances.csv`) | Sections 3.1, S5; Tables S12–S15; Figure S9 |

### `mopac_single_points/` – Sections 2.7 and 3.12, Figures 18 and 19

`B3mono2_sp/` (RM1_B3_h2_mono2) and `aqua_sp/` (RM1_CTRL_Eu_w9): RM1 single points on the optimized geometries (`sp.mop`, `sp.out`, `sp.arc`, `sp.aux`, and `sp.mgf` with the orbital coefficients). `code/orbital_analysis.py` reads them.

### `diagnostic_runs/` – exploratory runs cited in Sections 2.3–2.4

| File | Purpose |
|---|---|
| `diagnostic_runs_summary.csv` | one row per exploratory job (202 jobs of three early protocols): keywords, status, first-sphere composition, shortest Eu–O distance, proton-transfer counters |
| `diagnostic_runs_final_geometries.xyz` | last geometry of each of those jobs |
| `eu_pectin_campaign_v1.tar.gz`, `_v2`, `_v5` | input packages of the three protocols |

They support three statements of the Methods: COSMO gave artificially short Eu–O contacts and optimizer failures (v1); Sparkle/PM7 as distributed let hydroxyl and alkoxide-like oxygens collapse onto the metal (v2); water molecules transferred protons when the O–H lengths were not fixed (v5).

### `data/` – processed data behind the tables and figures

| File | Purpose | Used for |
|---|---|---|
| `results_v2_all_models.csv` | descriptor table, one row per optimization (388 rows, 79 columns; dictionary in Section 5) | Tables 3–6, S1, S2, S4, S5; Figures 11, 12, 15, 16, S5 |
| `final_geoms_v2_all.json` | final geometry of every optimization: elements, coordinates (Å, MOPAC atom order), heat of formation, net atomic charges | all structural figures; input of `reproduce_descriptors_v2.py` |
| `final_geometries_v2_all.xyz` | the same 388 structures as one multi-structure XYZ file | viewing in any molecular viewer |
| `manifest_v2_all.json` | construction record of the 388 optimizations (as `mopac_optimizations/manifest.json`, all models and the `ester2` starts) | residue atom map of the descriptors |
| `TableS9_integrity.csv` | integrity criterion applied to every optimization: broken and new heavy-atom bonds, loss of CO₂ (residue, C5···C6 distance, geometry and charge of the detached unit), net atomic charge of Eu | Section 2.5, Figure 6, Tables S2, S9 |
| `scans_semiempirical_v2.json` | rigid Eu–O scans of the four models: grid `r` (1.2–3.6 Å) and heats of formation of Eu³⁺, of each rigid ligand and of the complex | Figure 3, Tables 2, S3a–d |
| `check_fallback_vs_scans.csv` | point-by-point check of the Eu–O pair term (difference between two scans against the analytic difference of Eq. 1) | Section 2.3, Table S7 |
| `fig05b_calibration_map.csv` | position of the Eu–OH₂ and Eu–OCH₃⁻ scan minima on a grid of the two pair parameters | Figure 4b |
| `TableS8_R1_diagnostics.csv` | heats of formation of the 17 diagnostic single points | Table S8 |
| `phipsi_map.json` | MMFF94 φ/ψ maps of the disaccharide | Figure 2 |
| `TableS6_smiles.csv` | stereo-explicit SMILES as built, RDKit canonical SMILES, InChIKey, formula and charges of the ten hexasaccharides, the monomer controls and the scan ligands | Table S6, Section 2.2 |
| `TableS6_smiles_check_per_job.csv` | result of the comparison of each starting structure with its SMILES | Section 2.2 |
| `vibrational_summary_v2.csv` | constrained harmonic analysis of 229 structures (the 226 of the analysis set of the three production models and the three superseded RM1 M6 structures): number of constraints, gradient in the tangent space, imaginary frequencies, lowest frequencies, sum-rule violation of the Hessian | Section 3.11, Figure 17a, Table S10 |
| `displacement_test_results.csv` | displacement test, one row per re-optimization (186): status, energy change, RMSD from the parent structure, first-sphere changes and descriptors | Section 3.11, Figure 17b,c, Table S11 |
| `tables_v2/` | reference output of `reproduce_tables_v2.py`: `Table2_calibration.csv` … `Table6_matched_pairs.csv`, `TableS2_excluded.csv`, `TableS3_EuO_scans.csv`, `TableS4_controls.csv`, `TableS5_sensitivity.csv` and `key_numbers_v2.csv` (524 numbers quoted in the text, each with the section it belongs to) | Tables 2–6, S2–S5 and the text |
| `xray_benchmark/R10_COD_screening.csv` | the 1711 entries of the Crystallography Open Database that were screened, with the filters and the 17 selected | Section S5 |
| `xray_benchmark/cif/` | the 17 CIF files as downloaded from the COD (24 November 2026) | Table S12 |
| `xray_benchmark/xray_benchmark_results.csv`, `xray_benchmark_distances.csv` | benchmark results: one row per optimization; one row per Eu–O distance and model | Figure S9 |
| `xray_benchmark/TableS13_benchmark_by_model.csv`, `TableS14_benchmark_by_donor.csv`, `TableS15_benchmark_by_structure.csv` | Tables S13–S15 | Sections 3.1, S5 |
| `xray_benchmark/xray_benchmark_final_geometries.xyz` | optimized geometries of the 68 benchmark calculations | viewing |

### `code/`

| Script | Purpose | Article |
|---|---|---|
| `reproduce_descriptors_v2.py` | (A) compares status, final geometry, heat of formation, gradient and charge of Eu of every optimization output with the archived values; (B) recomputes all descriptors and the integrity table from the archived geometries; (C) checks the displacement test against its outputs | Tables S1, S9, S11 |
| `reproduce_tables_v2.py` | recomputes Tables 2–6, S2–S5 and the key numbers | Sections 2.5–3.11 |
| `analyze.py`, `sphere.py`, `residue_descriptors.py` | modules: run status from the `.out`, classification of every oxygen by connectivity, first-sphere descriptors (cutoff 3.0 Å), residue-level motifs | Section 2.6, S2 |
| `integrity_check.py` | module: bond criterion of the integrity test | Section 2.5 |
| `mopio.py` | module: MOPAC input writer (fixed O–H lengths) and output readers | Section 2.4 |
| `galA.py`, `assemble.py` | α-D-GalpA residue with checked stereochemistry (RDKit); assembly of the α-(1→4) hexasaccharides with prescribed φ/ψ and restrained MMFF94 relaxation | Section 2.2 |
| `place_eu.py` | placement of Eu and packing of the water molecules | Section 2.4 |
| `make_campaign.py` | generator of the 180 inputs of the original campaign | Sections 2.2–2.5 |
| `smiles_table.py` | writes Table S6 and compares every starting structure with its SMILES | Section 2.2, Table S6 |
| `phipsi_map.py` | MMFF94 φ/ψ map of the disaccharide | Figure 2 |
| `ligands.py`, `semi_scans.py` | rigid ligands and approach vectors; rigid scans with Sparkle/PM7, Sparkle/PM7+EuO and RM1 (runs MOPAC) | Section 2.3, Figure 3 |
| `check_scans.py` | compares the Sparkle/RM1 scans with the outputs of `revision/R4*` | Table S3c |
| `calibration_map.py` | check of the Eu–O pair term and map of the two-parameter calibration | Section 2.3, Figure 4, Table S7 |
| `constrained_modes.py` | harmonic analysis on the constrained surface from a MOPAC `FORCE` calculation; the command line analyses the outputs present and compares them with `data/vibrational_summary_v2.csv` | Sections 2.5, S4 |
| `xray_benchmark.py` | `build`: extracts the complexes from the CIF files and writes the 68 inputs; `analyze`: compares optimized and crystal Eu–O distances and writes Tables S13–S15 | Section S5 |
| `orbitals.py`, `orbital_analysis.py` | module that reads the basis and the orbital coefficients; HOMO/LUMO energies, Eu weights and net charge of Eu from the two single points | Section 3.12 |
| `summarize_diagnostic_runs.py` | wrote the two summary files of `diagnostic_runs/` from the raw exploratory outputs (kept for transparency; those outputs are not available) | Sections 2.3–2.4 |
| `run_campaign.ps1` | Windows PowerShell runner used for every campaign | Section 2.5 |

---

## 4. Reproduction

All commands are run from `code/`. Steps A–H need Python only and read the files of this package; steps I–K need MOPAC 23.2.5.

| Step | Command | Expected result |
|---|---|---|
| A. Outputs → descriptors | `python reproduce_descriptors_v2.py --compare` | outputs found for 388 of 388 optimizations; 360 converged, all equal to the archive; 186 of 186 displacement outputs; `identical values = True` |
| B. Tables and numbers | `python reproduce_tables_v2.py --out tables_check` | files equal to `data/tables_v2/` |
| C. SMILES | `python smiles_table.py --out smiles_check` | `start geometries checked: 172  matching: 172  not applicable: 8`; files equal to `data/TableS6_*.csv` |
| D. Pair term and calibration map | `python calibration_map.py --compare` | both tables equal to the archived ones (largest residual 0.045 kcal mol⁻¹ over 88 points) |
| E. Sparkle/RM1 scans | `python check_scans.py` | 139 single points, largest deviation 0.00000 kcal mol⁻¹ |
| F. Vibrational analysis | `python constrained_modes.py --runs ../revision/R3_force ../revision/R8_force_new_series --data ../data` | 9 FORCE jobs analysed, deviations 0 |
| G. Crystal benchmark | `python xray_benchmark.py analyze --runs ../revision/R10_xray_benchmark --out xray_check` ; `python xray_benchmark.py build --cif ../data/xray_benchmark/cif --out R10_rebuilt` | tables equal to `data/xray_benchmark/`; inputs equal to those of `revision/R10_xray_benchmark` |
| H. Frontier orbitals | `python orbital_analysis.py` | Eu weight in the LUMO 0.694 and 0.753; net charge of Eu −0.329 and +0.366 |
| I. Scans (MOPAC) | `python semi_scans.py` | heats of formation of three models; compare with `data/scans_semiempirical_v2.json` |
| J. Inputs (RDKit) | `python make_campaign.py <new_folder>` ; `python phipsi_map.py` | see limitation 4 |
| K. Any calculation (MOPAC) | `mopac <job>.mop` inside the job folder, or `powershell -ExecutionPolicy Bypass -File run_campaign.ps1 -Root <folder> -CampaignDir <subfolder> -MopacExe <path\to\mopac.exe>` | the runner skips job names ending in `_r<number>` and terminates running MOPAC processes when it starts |

**Dependencies.** Python ≥ 3.9 with numpy, pandas and scipy (steps A, B, D, E, F, H); RDKit (C, J and `build` of G); gemmi (`build` of G). Tested with Python 3.11.18, numpy 2.4.6, pandas 3.1.0, scipy 1.17.1, RDKit 2026.09.1, gemmi 0.7.5. MOPAC 23.2.5, open-source release (https://github.com/openmopac/mopac/releases); `semi_scans.py` calls the executable named in the environment variable `MOPAC_EXE` (default `mopac`).

---

## 5. Data dictionary of `results_v2_all_models.csv`

Distances in Å, energies in kcal mol⁻¹, angles in degrees; first sphere = oxygens within 3.0 Å of Eu.

| Column(s) | Definition |
|---|---|
| `job`, `model`, `series`, `pname`, `pattern`, `conf`, `mode`, `charge`, `sys` | identification; `charge` = net charge of the cluster; `sys` = job name without the model |
| `stage` | restart stage of the runner |
| `status` | `converged`: gradient norm below 1.0 kcal mol⁻¹ Å⁻¹; `cycles`: the 3000-cycle budget was exhausted; `ended`: self-consistence not achieved |
| `heat`, `gnorm`, `ncycles` | heat of formation (comparable only within identical composition), final gradient norm, geometry cycles |
| `CN`, `n_pectin`, `n_water`, `f_pectin`, `f_water` | coordination number; pectin and water oxygens; fractions of CN |
| `n_carboxylate`, `n_hydroxyl`, `n_ring_or_glycosidic`, `n_carbonyl_ester`, `n_ether_ester`, `n_acid_OH`, `n_carbonyl_acid`, `n_alkoxide`, `n_hydroxide`, `n_oxonium` | first-sphere oxygens by class |
| `d_pectin_mean`, `d_pectin_min`, `d_water_mean`, `d_water_min`, `d_water_max`, `d_carbox_mean`, `d_all_mean`, `d_all_std`, `d_min_any` | Eu–O distances by class; `d_min_any` = shortest Eu–O distance to any oxygen |
| `n_pectin_3.6`, `n_water_3.6` | pectin and water oxygens within 3.6 Å |
| `n_water_defects`, `n_acidOH`, `collapsed` | water-derived oxygens with a wrong number of hydrogens; acid O–H groups in the cluster; 1 if a first-sphere oxygen is closer than 2.0 Å |
| `n_residues_bound`, `residues_bound`, `n_carboxylates_bound`, `n_bidentate_carboxylate`, `n_monodentate_carboxylate`, `n_esters_bound`, `n_acids_bound`, `n_backbone_O`, `n_hydroxyl_O` | residue-level counts (index 0 = non-reducing end) |
| `motif` | residue-level motif, e.g. `A1:mono+O5 A2:mono` |
| `chain_rmsd`, `torsions_final`, `max_dtorsion`, `eu_min_chainO_start` | deformation of the chain relative to the starting structure; shortest Eu–pectin oxygen distance in the start |
| `intact`, `n_broken_bonds`, `n_new_bonds`, `n_CO2_lost` | integrity criterion: heavy-atom connectivity equal to that of the starting structure; number of broken and new bonds; carboxylates lost as CO₂ |
| `q_Eu` | net atomic charge of Eu |
| `superseded` | True for the six `ester` starts of M6 of the original campaign |
| `analysis_set` | converged, intact and not superseded (304 rows: 92 Sparkle/PM7+EuO, 99 Sparkle/RM1, 35 RM1, 78 Sparkle/PM7) |
| `vib_*` | constrained vibrational analysis: number of imaginary modes (all, below 20i, below 50i), two lowest frequencies, gradient in the tangent space, sum-rule violation, character of the lowest mode |
| `disp_*` | displacement test of the structure: re-optimizations converged, range of the energy change, largest RMSD and first-sphere Eu–O change, first sphere unchanged |

---

## 6. Provenance

**Molecular models and inputs.** `make_campaign.py` built the hexasaccharides from the ⁴C₁ residue of `galA.py`, assembled them with the φ/ψ of each conformer (`assemble.py`), placed Eu and the water molecules (`place_eu.py`) and wrote the inputs (`mopio.py`). The inputs of the other models (stages R2, R5) are the same starting geometries with the keywords of each model; those of R6 are new starts built on the same chains.

**Calculations.** All with MOPAC 23.2.5 for Windows, run by the author with `run_campaign.ps1`. Original campaign: 180 jobs; 23 of the 25 runs that failed to reach self-consistence were repeated (`_rerun`). A repetition replaces the original run unless it failed again. Totals over the four models: 388 optimizations, 360 converged, 20 exhausted the cycle budget, 8 did not reach self-consistence.

**Processed data.** Geometries, heats of formation and charges were read from the MOPAC outputs of each run; descriptors were computed from those geometries with the modules of `code/`. Step A of Section 4 repeats the comparison and the computation.

**Crystal structures.** Crystallography Open Database (https://www.crystallography.net/cod/; Gražulis et al., J. Appl. Crystallogr. 42 (2009) 726–729). The 17 structures were selected before any calculation, by the filters recorded in `R10_COD_screening.csv`; identifiers and the DOI of each structure report are in `R10_structures.csv` and in Table S12. Hydrogen atoms of the water molecules of one structure (XB03) were added.

**Exploratory runs.** Only the summary, the last geometries and the input packages are available (`diagnostic_runs/`).

---

## 7. Limitations and what is not in the package

1. **Auxiliary files (`.aux`) of the optimizations are not included** (they are most of the size of the raw data and can be regenerated by running the inputs). The package keeps, for every run, the input, the output and the archive file. The final geometry printed in the `.out` file has four decimals; step A confirms that it agrees with the archived geometry within 5 × 10⁻⁵ Å for the 360 converged runs.
2. **Hessian outputs are included for nine representative structures only**: [Eu(H₂O)₉]³⁺, the galacturonate chelate and a B3 hexasaccharide complex, each with the three production models. The inputs of all `FORCE` calculations are included and reproduce the others; the results of all of them are in `data/vibrational_summary_v2.csv`.
3. **Runs that did not converge.** For the 28 runs that did not reach the convergence criterion the last geometry is not printed in the `.out` file, so their archived geometries cannot be compared with the outputs here. They are excluded from every analysis and listed in Table S2.
4. **Starting geometries are not bit-reproducible by re-running the construction code.** With a different RDKit version `make_campaign.py` gives the same atoms in the same order but slightly different coordinates (MMFF relaxation, placement of Eu, packing of the water molecules), and the fully ionized φ/ψ map shifts by one grid point. The archived inputs are the reference.
5. **Scans.** The outputs of the scans of Sparkle/PM7, Sparkle/PM7+EuO and RM1 were not kept; `semi_scans.py` regenerates them (identical for water, methanol, methoxide and acetate; within 0.001 kcal mol⁻¹ for ether and ester).
6. **Benchmark.** The `analyze` step reads geometries with four decimals; errors are reproducible to 10⁻⁴ Å.
7. **Energies.** Heats of formation are comparable only within series of identical composition.
8. **Not included:** the scripts that draw the figures, the scripts that generated the inputs of the revision stages from the outputs of earlier stages (the inputs themselves are included), intermediate files, and superseded versions of data and documents. They are available from the author on request.
