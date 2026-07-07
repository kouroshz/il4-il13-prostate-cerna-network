# Quentin Verification TODO

## Priority 1: CircRNA Enrichment Panels

Verify or reproduce the Figure 2 circRNA host-gene enrichment panels using the exact input/background files documented in `documentation/final/quentin_circrna_enrichment_handoff.md`.

Key question: should the manuscript retain, revise, move, or remove the existing circRNA enrichment panels given that the current reproducible gprofiler2 run found zero adjusted-significant circRNA host-gene terms?

## Priority 2: Figure 4D Network Source

Confirm that `IL4 IL13 Integrated ceRNA Network.xlsx` is the curated source for Figure 4D and the pro-fibrotic CEMIP/TNC/LIFR-filtered network.

Current public TSV copies:

- `data/processed/final/integrated_ceRNA_network_curated.tsv`
- `results/final/source_data/Supplementary_Table_integrated_ceRNA_network.tsv`
- `results/final/networks/integrated_ceRNA_network_curated.tsv`

## Priority 3: Figure 4A/B/C Source

Confirm whether Figure 4A, Figure 4B, and Figure 4C also came from the integrated workflow, related treatment-specific tables, or a separate Cytoscape/export workflow.

## Priority 4: Candidate-Axis Score And Energy Values

Confirm score and predicted binding-energy values for:

- miR-140-3p -> CEMIP: score 152, energy -16.17
- miR-140-3p -> circPAPPA: score 147, energy -18.27
- miR-135b-5p -> TNC: score 152, energy -11.55
- miR-135b-5p -> circLIFR: score 143, energy -9.43
- miR-625-3p -> LIFR: score 157, energy -14.26

## Priority 5: Coordinate Interpretation

The workbook has columns named `mir.strt`, `mir.stp`, `gene.strt`, and `gene.stp`. Candidate `mir.strt` values exceed mature miRNA lengths, so current public wording does not call these coordinate fields.

Please confirm whether these columns represent:

- target coordinates,
- miRNA coordinates,
- row/index plus coordinate metadata,
- miRanda output fields with shifted interpretation,
- malformed parse,
- or unknown.

Until confirmed, public wording should report miRanda scores and predicted binding energies only.

## Priority 6: Candidate Labels

Confirm final labels for:

- miR-140-3p / CEMIP / circPAPPA
- miR-135b-5p / TNC / circLIFR
- miR-625-3p / LIFR

## Priority 7: Experimental-Result Status

Confirm the final experimental-result status and panel references for the three candidate axes:

- miR-140-3p/CEMIP: strongest functional support; pursued further.
- miR-135b-5p/TNC/circLIFR: evaluated candidate; negative or less supportive result unless updated by authors.
- miR-625-3p/LIFR: evaluated candidate; negative or less supportive result unless updated by authors.
- circPAPPA: predicted circRNA partner in the miR-140-3p/CEMIP axis; not functionally supported if circPAPPA knockdown remains negative.
