# Compleasm fork for Binny2

Upstream base: `9d83c884ea977b52d99aab24417b6555f2ea5dd3` (Compleasm 0.2.9, `#68`).

Fork version: `0.2.9+binny2.1` in `_version.py`. Pin the fork commit, not only this version string. `__upstream_base__` records the upstream commit above.

## Candidate-loss bug

Genome-mode `load_hmmsearch_output` set `best_protein = None` for each HMMER domtblout file. After the first reference protein name, `if best_protein is not None and best_protein != protein_name: continue` dropped every later reference protein. The stored name was the first name seen. A failing first protein still occupied that slot, so a later protein that passed the BUSCO score cutoff was discarded. The comment described a best-score filter. The code implemented a first-seen filter.

Same reference protein at two genomic locations was already retained, because the filter keyed only `protein_name`. The loss is alternate reference proteins, which can sit on different physical loci. Protein mode parses domtblout separately and already kept every passing target. It is unchanged.

## Patch

`load_hmmsearch_output` groups rows by the full target `protein_name|contig:start-end`, then by HMM domains. A candidate is kept when its full-sequence bitscore (domtblout column 8) is greater than or equal to the BUSCO cutoff. Domain rows of one target are merged with the existing to-from interval arithmetic (`hmm_to - hmm_from`, no inclusive `+1`). Overlapping domains extend one interval. They do not become extra loci.

Returned mapping keys stay `protein_name|contig:start-end`. `MiniprotAlignmentParser.Run` still classifies with `Ost_eval`. It also writes `hmm_candidate_evidence.tsv` beside `full_table.tsv`.

HMM filenames are sorted. Evidence rows and BUSCO ids written to the full table are sorted. Within one BUSCO, rows are ordered by `I+L` descending, then protein name and genomic coordinates (`mergesort`). Equal `I+L` no longer follows Miniprot row order when choosing the classified representative.

Miniprot is still invoked with the default `--outs=0.95`. This fork does not change that threshold.

## Sidecar schema 1

Path: `<output>/<lineage>/hmm_candidate_evidence.tsv`.

One row is one HMM-supported reference protein at one genomic locus. Columns: `schema_version`, `busco_id`, `protein_name`, `locus_id`, `contig`, `genomic_start`, `genomic_end`, `strand`, `hmm_score`, `hmm_qlen`, `hmm_tlen`, `hmm_matched_length`, `hmm_domains`, `hmm_domain_count`, `miniprot_rank`, `miniprot_identity`, `miniprot_positive`, `miniprot_score`, `source_hmm_file`.

`locus_id` is `contig:start-end` from the translated-protein header. Contig names may contain colons; coordinates are the trailing field. `hmm_domains` is `hmm_from-hmm_to`, sorted. Header comments record schema version, fork version, and upstream base.

Two different reference proteins at the same coordinates remain two rows. Binny2 can collapse those to one physical locus. Two HMM-supported mappings at different coordinates remain two rows.

## Classification

`full_table.tsv` is still Compleasm's classification, not the raw locus list. Different reference proteins that are each a single complete hit can still be labeled `Single`, and the table keeps the higher `I+L` protein. The sidecar keeps both observations. A previously dropped protein with a higher `I+L` can change which reference is reported, and a previously dropped close alternate can change `Single` to `Duplicated`. Those differences are corrections of the first-candidate filter. They are not forced by collapsing loci inside Compleasm.

## Tests

`python3 -m unittest tests.test_hmm_candidate_retention -v`

The repository has no upstream unit suite. The new tests cover two passing candidates, row-order permutations, a failing candidate placed first, multi-domain aggregation, same locus with two references, distinct loci, three passing candidates, directory order, and sidecar row order.
