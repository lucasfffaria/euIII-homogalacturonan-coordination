# Supplementary Material

**Semiempirical exploration of Eu(III) coordination in microsolvated homogalacturonan fragments: effects of the methylesterification pattern**

This folder is the audit and reproducibility package of the article. Together with the manuscript, it lets an independent reader:

1. inspect every MOPAC input and output behind the 180 optimizations (`mopac_optimizations/`);
2. rebuild the descriptor table (Table S1) and the final geometries from those raw outputs (`code/reproduce_descriptors.py`);
3. recompute Tables 2–6, Tables S2–S5 and the quantitative statements of the text from the descriptor table (`code/reproduce_tables.py`);
4. repeat the Eu–O calibration scans, the RM1 single points of Section 3.10 and the construction of the inputs (`code/`);
5. examine the exploratory runs on which three statements of the Methods rest (`diagnostic_runs/`).

Every file has a single location. The only repeated contents are the per-job files that MOPAC needs in each job folder (`eu_fix.txt`).

---

## 1. Quick check (Python only, about 2 minutes)

```bash
cd code
python reproduce_descriptors.py --compare ../data   # raw MOPAC outputs -> Table S1 + final geometries (~1 min, reads ~1.7 GB)
python reproduce_tables.py --data ../data           # Table S1, scans, geometries -> Tables 2-6, S2-S5, key numbers (seconds)
```

The first command must end with `identical bytes = True` for `results_final_v2.csv`, `final_geoms_final.json` and `final_geometries_all.xyz`. The second writes one CSV per table to `./tables_out` (Section 4).

---

## 2. Contents

### Documents

| File | Content |
|---|---|
| `Supplementary_Information.docx`, `.pdf` | Supplementary Information: S1 computational details, S2 descriptor definitions, Tables S1–S5, Figures S1–S3 |
| `README.md` | this file |

### `data/` – processed data (inputs of every table and figure)

| File | Content | Used for |
|---|---|---|
| `results_final_v2.csv` | descriptor table, one row per optimization (180 rows × 58 columns; dictionary in Section 3) | Table S1; source of Tables 3–6, S2, S4, S5 and Figures 4–10, S1–S3 |
| `final_geoms_final.json` | final structure of each optimization: `el` (elements), `xyz` (Å, MOPAC atom order), `charges` (Mulliken), `heat` (heat of formation, kcal mol⁻¹), `gnorm`, `ncycles` (cycles of the selected run), `heat_traj` (heat of formation along the optimization: every 5th point for original runs, every point for reruns); `stage`/`status` are present only in entries taken from original runs | Figures 4, 5, 8, S1, S2; cluster sizes (Section 2.5); Mulliken charges (Section 3.10) |
| `final_geometries_all.xyz` | the same 180 structures as one multi-structure XYZ file (comment line: job, model, pattern, conformer, start, status, heat) | viewing in any molecular viewer |
| `scans_semiempirical.json` | rigid Eu–O scans: distance grid `r` (1.2–3.6 Å) and, for each Hamiltonian (`Sparkle/PM7`, `Sparkle/PM7*`, `RM1`), heats of formation of Eu³⁺ (`Eu3+`), of each rigid ligand (`<ligand>_ligand`) and of the complex along `r` (`<ligand>`); binding energy = H(complex) − H(Eu³⁺) − H(ligand) | Table 2, Table S3, Figure 3 |
| `phipsi_map.json` | MMFF94 (ε = 4r) φ/ψ maps of the α-(1→4) disaccharide with both residues as free carboxylate (`A`), methyl ester (`M`) or acid (`H`): 15° grid, energies relative to the map minimum (kcal mol⁻¹), local minima below 6 kcal mol⁻¹ | Figure 2 (map `A`), Section 2.2 |

### `mopac_optimizations/` – raw data of the 180 optimizations (~1.7 GB)

One folder per job (names in Section 3). Each folder is self-contained, so the job can be re-run in place.

| File | Content |
|---|---|
| `<job>.mop` | MOPAC input: keyword line, title, geometry in mixed Cartesian/internal coordinates (O–H bonds as frozen internal coordinates) |
| `<job>_start.xyz` | starting geometry in MOPAC atom order (reference for `chain_rmsd`, `max_dtorsion`, `eu_min_chainO_start`) |
| `<job>.out`, `<job>.aux`, `<job>.arc` | MOPAC output; auxiliary output (geometry, charges, heat-of-formation trajectory); final archive (`.arc`, present only for runs that converged) |
| `eu_fix.txt` | the Sparkle/PM7* correction, read through `EXTERNAL=eu_fix.txt` (two lines: `XFAC_Eu O 2.207500`, `ALPB_Eu O 2.500000`); identical in the 101 Sparkle/PM7* folders, because MOPAC reads it from the job folder |
| `<job>_rerun.*` | repetition of a run that failed to reach self-consistence, from the same starting geometry with `SHIFT=15 PULAY CAMP ITRY=999` added (23 jobs; Section 5) |
| `manifest.json` | construction record of every job: model, charge, atom count, start mode, target residues and Eu–O target distance, waters around Eu and on carboxylates, chain torsions and MMFF energy, residue atom map in MOPAC atom order (used by the analysis), permutation from builder to MOPAC atom order |
| `progress_campaign6.log` | log of the production run (7 Sep 2026: 180 jobs in 158 min) |

### `mopac_single_points/` – Sections 2.7 and 3.10, Figure 11

`B3mono2_sp/` (RM1_B3_h2_mono2, q = 0) and `aqua_sp/` (RM1_CTRL_Eu_w9, [Eu(H₂O)₉]³⁺, q = +3): RM1 single points (`1SCF GRAPHF AUX VECTORS ALLVEC MULLIK BONDS`) on the final geometries, with `sp.mop`, `sp.out`, `sp.aux`, `sp.mgf` (orbital coefficients used for the grids) and `sp.arc`. These files were regenerated on 13 Oct 2026 with `code/fig_orbitals.py` and MOPAC 23.2.5 (Linux build), because the original single-point outputs had not been kept; they reproduce the panels of Figure 11 pixel for pixel and the numbers of Section 3.10.

### `diagnostic_runs/` – exploratory runs cited in Sections 2.3–2.4

| File | Content |
|---|---|
| `diagnostic_runs_summary.csv` | one row per job that produced output (202 jobs: v1 44, v2 88, v5 70): protocol, keyword line actually run, status, gradient, first-sphere composition, shortest Eu–O distance (`d_min_any`), `collapsed` (an O within 2.0 Å of Eu), proton-transfer counters (`n_hydroxide`, `n_oxonium`, `n_water_defects`, `n_alkoxide`, `n_acidOH`), motif |
| `diagnostic_runs_final_geometries.xyz` | last geometry of each of those jobs |
| `eu_pectin_campaign_v1.tar.gz`, `_v2`, `_v5` | original input packages (in v5, `LET` was added afterwards to 146 inputs that had not run yet; the summary records the keywords actually run) |

They document three statements of the Methods: v1 (Sparkle/PM7 as distributed, COSMO ε = 78.4, Cartesian) – COSMO gave artificially short Eu–O contacts and optimizer failures; v2 (Sparkle/PM7 as distributed, gas phase, O–H frozen) – hydroxyl and alkoxide-like oxygens collapsed onto the metal at 1.3–1.45 Å; v5 (Sparkle/PM7*, Cartesian, no constraint) – water transferred protons when unconstrained. The raw output folders of these campaigns (~1.7 GB, mostly `.aux` trajectories) were not retained; the summary and geometries were written from them by `code/summarize_diagnostic_runs.py`.

### `code/` – analysis, reproduction and input-generation code

| Script | Role | Article |
|---|---|---|
| `reproduce_descriptors.py` | rebuilds the three files of `data/` that derive from the optimizations (analysis of every run, merge of reruns, final geometries, residue descriptors); `--compare` checks them against `data/` | Table S1 |
| `reproduce_tables.py` | recomputes Tables 2–6, S2–S5 and `key_numbers.csv` (counts, ranges, correlations and differences quoted in the text) | Sections 2.5–3.10 |
| `analyze.py` | first-sphere descriptors of one run (cutoff 3.0 Å) and run status from the `.out` | Section 2.6 |
| `sphere.py` | classification of every oxygen by connectivity | Section 2.6, S2 |
| `residue_descriptors.py` | residue-level descriptors with the strict carboxylate definition (kind A only) | Section 2.6, S2 |
| `mopio.py` | MOPAC input writer (mixed internal/Cartesian, frozen O–H) and `.aux`/`.arc` readers | Section 2.4, S1 iii |
| `galA.py`, `assemble.py` | α-D-GalpA residue with CIP-checked stereochemistry; assembly of α-(1→4) oligomers with prescribed φ/ψ and restrained MMFF94 relaxation (RDKit) | Section 2.2 |
| `phipsi_map.py` | MMFF94 φ/ψ map of the disaccharide | Section 2.2, Figure 2 |
| `place_eu.py` | Eu placement (penalty minimization on the target Eu–O distances) and packing of the water molecules | Section 2.4 |
| `make_campaign.py` | generator of the 180 inputs (main series, hydration series, controls) | Sections 2.2–2.5, Table 1 |
| `run_campaign.ps1` | Windows PowerShell runner used for the production run and the reruns | Section 2.5 |
| `ligands.py`, `semi_scans.py` | rigid ligands and Eu approach vectors; MOPAC scans with the three Hamiltonians | Section 2.3, Table 2, Table S3, Figure 3 |
| `fig_orbitals.py`, `orbitals.py`, `render.py` | RM1 single point; orbitals and density on a grid from `.aux`/`.mgf` (Slater functions); Coulomb potential of the Mulliken charges; HOMO/LUMO energies, Eu weights and charge; PyVista panels | Sections 2.7, 3.10, Figure 11 |
| `summarize_diagnostic_runs.py` | wrote `diagnostic_runs/` from the raw exploratory outputs (kept for transparency; those outputs are no longer available) | Sections 2.3–2.4 |

### Where each result comes from

| Article item | Data | Code |
|---|---|---|
| Table 1, Sections 2.2 and 2.4 | `mopac_optimizations/manifest.json`, `*.mop`, `*_start.xyz` | `make_campaign.py`, `galA.py`, `assemble.py`, `place_eu.py`, `mopio.py` |
| Figure 2 | `data/phipsi_map.json` | `phipsi_map.py` |
| Section 2.3, Table 2 (scans), Figure 3, Table S3 | `data/scans_semiempirical.json` | `ligands.py`, `semi_scans.py`, `reproduce_tables.py` |
| Sections 2.3–2.4 (exploratory runs) | `diagnostic_runs/` | `summarize_diagnostic_runs.py` |
| Section 2.5 (protocol and counts), Table S2 | `mopac_optimizations/`, `data/results_final_v2.csv` | `run_campaign.ps1`, `reproduce_descriptors.py`, `reproduce_tables.py` |
| Section 2.6, Table S1 | `data/results_final_v2.csv` | `analyze.py`, `sphere.py`, `residue_descriptors.py`, `reproduce_descriptors.py` |
| Section 3.1, Table 2 (controls), Table S4 | `results_final_v2.csv` (`series = control`), `final_geoms_final.json` | `reproduce_tables.py` |
| Sections 3.2–3.8, Tables 3–6, Figures 4–10, S1–S3 | `results_final_v2.csv`, `final_geoms_final.json` | `reproduce_tables.py` |
| Section 3.9, Table S5 | `results_final_v2.csv` (`series = sens_hydration`) | `reproduce_tables.py` |
| Sections 2.7 and 3.10, Figure 11 | `mopac_single_points/`, `final_geoms_final.json` | `fig_orbitals.py`, `orbitals.py`, `render.py` |

---

## 3. Names and data dictionary

Job names are `<model>_<pattern>_<conformer>_<start>[_cw0|_cw4]` (e.g. `RM1_B3_h2_mono2`); controls are `<model>_CTRL_<system>_<variant>` (e.g. `PM7_CTRL_Eu_w9` = [Eu(H₂O)₉]³⁺).

| Code | Meaning |
|---|---|
| `PM7` | Sparkle/PM7* = Sparkle/PM7 with the Eu–O correction of `eu_fix.txt` ("PM7*" in the tables) |
| `RM1` | RM1 lanthanide model (Eu with a 5d/6s/6p basis) |
| pattern | codes of Table 1: A6, H3blk, H3alt, P3, B3, L3, B2, D2, S1, M6; column `pattern` gives the sequence (A = –COO⁻, M = –COOCH₃, H = –COOH, non-reducing → reducing end). In the figures H3blk and H3alt are labelled "A3H3 block" and "A3H3 alt." |
| conformer | `h3` 3₁-like helix (φ, ψ = 75°, 105°), `h2` 2₁-like helix (105°, 135°), `cl` irregular coil; `-` for controls |
| start (column `mode`) | `bid` bidentate, `mono2` monodentate plus bridge, `hyd` hydrated ([Eu(H₂O)₉]³⁺ away from the chain), `ester` ester carbonyl (M6); controls: `w8`–`w12` (number of waters), `bid`, `hyd`, `ester` |
| `_cw0`, `_cw4` | hydration-sensitivity series (Section 3.9, Table S5): 0 or 4 hydrogen-bonded waters on every free carboxylate not bound to Eu (main series: 2) |
| `series` | `main` (150 jobs), `sens_hydration` (12), `control` (18) |

Columns of `results_final_v2.csv` (distances in Å, energies in kcal mol⁻¹, angles in degrees; first sphere = oxygens within 3.0 Å of Eu):

| Column(s) | Definition |
|---|---|
| `job`, `model`, `series`, `pname`, `pattern`, `conf`, `mode`, `charge` | identification; `charge` = net charge q of the cluster |
| `stage` | restart stage of the runner (0 for every job) |
| `status` | `converged`: gradient norm below 1.0 kcal mol⁻¹ Å⁻¹; `cycles`: the 3000-cycle budget was exhausted ("cycles" in Table S1, "cycle budget exhausted" in Table S2); `ended`: self-consistence not achieved ("SCF") |
| `heat` | heat of formation of the analysed geometry (comparable only within identical composition) |
| `gnorm` | final gradient norm (kcal mol⁻¹ Å⁻¹) |
| `ncycles` | geometry cycles of the selected run (points of its heat-of-formation trajectory) |
| `CN`, `n_pectin`, `n_water`, `f_pectin`, `f_water` | coordination number; pectin and water oxygens; f_pectin = n_pectin/CN, f_water = n_water/CN |
| `n_carboxylate`, `n_hydroxyl`, `n_ring_or_glycosidic`, `n_carbonyl_ester`, `n_ether_ester`, `n_acid_OH`, `n_carbonyl_acid`, `n_alkoxide`, `n_hydroxide`, `n_oxonium` | first-sphere oxygens by class (`sphere.py`; bond criterion 1.2 × sum of covalent radii) |
| `d_pectin_mean`, `d_pectin_min`, `d_water_mean`, `d_water_min`, `d_water_max`, `d_carbox_mean`, `d_all_mean`, `d_all_std` | first-sphere Eu–O distances by class |
| `d_min_any` | shortest Eu–O distance to any oxygen |
| `n_pectin_3.6`, `n_water_3.6` | pectin / water oxygens within 3.6 Å |
| `n_water_defects` | water-derived oxygens with a wrong number of hydrogens (hydroxide, oxonium, oxide) |
| `n_acidOH` | acid O–H groups in the whole cluster (3 in H3blk/H3alt, 0 elsewhere unless a proton moved) |
| `collapsed` | 1 if a first-sphere oxygen is closer than 2.0 Å |
| `n_residues_bound`, `residues_bound` | residues with any oxygen within 3.0 Å (index 0 = non-reducing end) |
| `n_carboxylates_bound`, `n_bidentate_carboxylate`, `n_monodentate_carboxylate` | free carboxylates (kind A only) bound through one or both O6 |
| `n_esters_bound`, `n_acids_bound` | ester (M) or acid (H) residues bound through their C6 oxygens |
| `n_backbone_O`, `n_hydroxyl_O` | bound ring/glycosidic (O1, O4, O5) and hydroxyl (O2, O3) oxygens |
| `motif` | residue-level motif, e.g. `A1:mono+O5 A2:mono` (residue 1, kind A, monodentate carboxylate plus its O5; residue 2 monodentate); `bi`/`mono` carboxylate, `est`/`acid` ester or acid O6, `+On` other oxygens |
| `chain_rmsd` | heavy-atom RMSD of the chain relative to `<job>_start.xyz` (Kabsch superposition) |
| `torsions_final` | final φ/ψ of the five linkages (`phi/psi;...`) |
| `max_dtorsion` | largest change of any glycosidic torsion relative to the start |
| `eu_min_chainO_start` | shortest Eu–pectin oxygen distance in the start geometry |
| `rerun` | 1 if the row comes from the repetition of the run (`<job>_rerun.*`) |
| `scf_fail` | constant 0 (column kept from the original pipeline; SCF failures have `status = ended`) |

---

## 4. Reproduction

### Requirements

| Step | Needs | Tested with |
|---|---|---|
| A, B | Python ≥ 3.9, numpy, pandas, scipy | Python 3.10.12, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3 (A also with Python 3.13, pandas 3.0) |
| C | MOPAC 23.2.5, RDKit (ether and ester geometries) | MOPAC 23.2.5 Linux, RDKit 2026.09.1 |
| D | MOPAC 23.2.5, numpy, scipy, PyVista | MOPAC 23.2.5 Linux, PyVista 0.49.1 |
| E, F | RDKit, scipy | RDKit 2026.09.1 |
| G | MOPAC 23.2.5 (Windows PowerShell for the runner) | MOPAC 23.2.5 Windows, 20-core Intel i7-14700 |

MOPAC 23.2.5 is the open-source release (https://github.com/openmopac/mopac/releases). `semi_scans.py` and `fig_orbitals.py` call the executable named in the environment variable `MOPAC_EXE` (default: `mopac` on the PATH); `run_campaign.ps1` takes `-MopacExe <path>`. All commands below are run from `code/`.

**A. Tables and numbers of the text (seconds).**
`python reproduce_tables.py --data ../data --out ../tables_out` writes `Table2_calibration.csv`, `Table3_by_start.csv`, `Table4_by_pattern.csv`, `Table5_waters.csv`, `Table6_matched_pairs.csv`, `TableS2_not_converged.csv`, `TableS3_EuO_scans.csv`, `TableS4_controls.csv`, `TableS5_sensitivity.csv` and `key_numbers.csv`. Definitions (Sections 2.5–2.6 and S2): converged = `status == "converged"`; main series = `series == "main"`; inner-sphere starts = `bid`, `mono2`, `ester`; exclusion rule CN < 6 (one structure, PM7_M6_cl_ester).

**B. Descriptor table and final geometries from the raw outputs (~1 min).**
`python reproduce_descriptors.py --runs ../mopac_optimizations --out ./descriptors_out --compare ../data` also writes the intermediate tables `results_original_runs.csv` (original runs only) and `results_merged.csv` (after the merge, before the residue descriptors).

**C. Eu–O scans (a few minutes).**
`python semi_scans.py` writes the inputs/outputs to `./scan_runs` and the heats of formation to `./scans_semiempirical.json`; compare with `data/scans_semiempirical.json`.

**D. Single points and Figure 11 panels.**
`python fig_orbitals.py ../data/final_geoms_final.json:RM1_B3_h2_mono2 0 B3mono2` and `python fig_orbitals.py ../data/final_geoms_final.json:RM1_CTRL_Eu_w9 3 aqua` run the single points in `B3mono2_sp/` and `aqua_sp/`, print the HOMO/LUMO energies, the Eu weight of HOMO and LUMO and the Eu Mulliken charge, and render the panels (`<prefix>_homo.png`, `_lumo.png`, `_density_esp.png`).

**E. φ/ψ map.** `python phipsi_map.py` writes `phipsi_map.json` (a few minutes; see limitation 1).

**F. Optimization inputs.** `python make_campaign.py <new_folder>` writes the 180 job folders and `manifest.json` (see limitation 1).

**G. Optimizations.** Each job runs in its folder with `mopac <job>.mop` (Sparkle/PM7* jobs need `eu_fix.txt` in the same folder). The production run used `powershell -ExecutionPolicy Bypass -File run_campaign.ps1 -Root <folder> -CampaignDir <subfolder> -MopacExe <path\to\mopac.exe>` with 13 jobs in parallel and `THREADS=2` (158 min for the 180 jobs). The runner terminates every running MOPAC process when it starts (see the NOTE in the script).

---

## 5. Provenance

**Construction (Sections 2.2–2.4).** `make_campaign.py` built the hexasaccharides from the CIP-verified ⁴C₁ residue (`galA.py`), assembled them with the φ/ψ of each conformer and relaxed them with restrained MMFF94 (`assemble.py`), then placed Eu and the water molecules (`place_eu.py`) and wrote the inputs (`mopio.py`). As encoded in the inputs and in `manifest.json`:

- target Eu–O distance of the starting contacts 2.48 Å;
- nine waters around Eu: on the free vertices of a tricapped trigonal prism at 2.48 Å (at most nine minus the number of starting chain contacts), the others in a second shell at 4.6 Å, or farther when there is no room (`n_first_shell_waters` gives how many are in the first shell);
- main series and galacturonate controls: in addition, two hydrogen-bonded waters on each free carboxylate not bound to Eu at the start (`n_carboxylate_waters`), 9–21 waters in total (`n_water`); hydration series: 0 (`_cw0`) or 4 (`_cw4`) such waters per free carboxylate, for B3 and L3 with both Hamiltonians and for A6 and P3 with Sparkle/PM7*;
- every O–H bond is an internal coordinate with frozen length: 0.96 Å for water; the sugar hydroxyl and acid O–H keep their MMFF94 lengths (0.97–0.98 Å);
- cluster sizes: 145–205 atoms (hexasaccharide clusters) and 25–56 atoms (controls);
- RM1 was not run for the net-anionic A6 and P3.

The inputs were first written for an earlier campaign (v4, `CYCLES=8000`) that was stopped; the production inputs have the same coordinates and start files, with the production keywords:
`PM7 CHARGE=q AUX GNORM=1.0 CYCLES=3000 THREADS=2 GEO-OK T=3D DUMP=30M BONDS EF LET EXTERNAL=eu_fix.txt` (Sparkle/PM7*) and the same line with `RM1` and without `EXTERNAL` (RM1).

**Production run (Section 2.5).** 180 jobs, 7 Sep 2026, 13 in parallel (`progress_campaign6.log`). The runner's automatic restart was not triggered (`stage` = 0 for every job). Outcome of the original runs: 145 converged (90 of 101 Sparkle/PM7*, 55 of 79 RM1); 10 Sparkle/PM7* runs used the whole 3000-cycle budget (residual gradient norms 1.5–9.5 kcal mol⁻¹ Å⁻¹); 25 failed SCF (24 RM1 and PM7_A6_cl_hyd).

**Repetitions.** 23 of the 25 runs that failed to reach self-consistence (22 RM1 runs and PM7_A6_cl_hyd) were repeated from their starting geometries with `SHIFT=15 PULAY CAMP ITRY=999` added and 3000 cycles (`<job>_rerun.*`): 18 of the RM1 repetitions converged, 4 failed again, and the repetition of PM7_A6_cl_hyd exhausted the cycle budget (final gradient 6.8). RM1_D2_cl_bid and RM1_D2_cl_mono2 have no repetition, and the cycle-limited Sparkle/PM7* runs were not continued.

**Merge (`reproduce_descriptors.py`).** A repetition replaces the original run unless it failed again (status `ended`). The table therefore uses 19 repetitions (18 converged RM1 runs and PM7_A6_cl_hyd) and the original run of the other 161 jobs: 163 converged structures (90 Sparkle/PM7*, 73 RM1) and 17 that did not reach the criterion and are excluded from all statistics (descriptors of their last geometry) – 11 with `status = cycles` (10 Sparkle/PM7* runs and the repetition of PM7_A6_cl_hyd) and 6 RM1 runs with `status = ended`.

**Descriptors (Section 2.6).** `analyze.py` (oxygen classes, first sphere < 3.0 Å, distances, chain RMSD and torsions against the start) and `residue_descriptors.py` (strict carboxylate counting) applied to the final `.aux` geometry of each selected run.

**Scans (Section 2.3).** `semi_scans.py` with the rigid ligands of `ligands.py` (water, methanol, methoxide and bidentate acetate with fixed geometries; dimethyl ether and methyl acetate MMFF-optimized with RDKit): 1SCF calculations along `r` with Sparkle/PM7, Sparkle/PM7* and RM1. The original scan outputs were not kept; the re-run of 13 Oct 2026 reproduces `data/scans_semiempirical.json` exactly for water, methanol, methoxide and acetate and within 0.001 kcal mol⁻¹ for ether and ester.

**Single points (Section 2.7).** See `mopac_single_points/`.

---

## 6. Verification performed (13–14 Oct 2026)

- `reproduce_descriptors.py`, run on the folders of this package, rebuilds `results_final_v2.csv`, `final_geoms_final.json` and `final_geometries_all.xyz` byte for byte.
- `reproduce_tables.py` reproduces the 306 values of Tables 2–6 and of the text checked against the article and all 1885 cells of Tables S1–S5 of the Supplementary Information; its output files are byte-identical under Python 3.10/pandas 2.3 and Python 3.13/pandas 3.0.
- The scans and single points regenerated with MOPAC 23.2.5 agree with the archived data (Section 5); the Figure 11 panels rendered from the regenerated single points are pixel-identical to the published ones.
- All figures of the article (1–11) and of the SI (S1–S3) were regenerated from the files in `data/` (and `mopac_single_points/` for Figure 11) and match the published versions.

---

## 7. Limitations

1. **Start geometries are not bit-reproducible by re-running the construction code.** With RDKit 2026.09.1, `make_campaign.py` produces the same atoms in the same order but different coordinates (MMFF relaxation of the charged residues, Eu placement, random packing of the waters). The archived inputs in `mopac_optimizations/` are the reference. Likewise, `phipsi_map.py` reproduces the minima of the `M` and `H` maps, but the `A` map (fully ionized) shifts with the RDKit version (global minimum at 105°/150° instead of 105°/135°).
2. **Regenerated outputs.** The original outputs of the scans and of the two single points were not kept; those in this package were regenerated as described in Section 5.
3. **Exploratory runs.** Only the summary, the last geometries and the input packages of the exploratory campaigns are available.
4. **Energies.** Heats of formation are comparable only within series of identical composition (Section 2.3 of the article).
5. **Not included:** the plotting scripts used to draw the figures and the scripts that assembled the manuscript and the SI documents; the figures can be redrawn from `data/` with any plotting tool, and the scripts are available from the authors on request.
