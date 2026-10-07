"""Regression tests for HMM-supported BUSCO candidate retention.

Upstream ``load_hmmsearch_output`` keeps the first acceptable ``protein_name`` in
each domtblout file and drops every later reference protein, so retained loci
depend on HMM row order. These tests require every passing candidate to survive,
independent of row order and directory order.
"""

import os
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import compleasm  # noqa: E402


BUSCO = "10at2"
QUERY = BUSCO
CUTOFF = {BUSCO: 40.0}


def _domtbl_line(target, score, hmm_from, hmm_to, query=QUERY, tlen=300, qlen=200):
    columns = ["-"] * 22
    columns[0] = target
    columns[2] = str(tlen)
    columns[3] = query
    columns[5] = str(qlen)
    columns[7] = "{:.2f}".format(score)
    columns[15] = str(hmm_from)
    columns[16] = str(hmm_to)
    return " ".join(columns)


def _target(reference, contig, start, stop):
    protein = "{}_{}_0:001".format(BUSCO, reference)
    location = "{}:{}-{}".format(contig, start, stop)
    return "{}|{}".format(protein, location), protein, location


def _write_folder(files):
    """files: list of (filename, list of lines)."""
    directory = tempfile.mkdtemp(prefix="compleasm-hmm-")
    for name, lines in files:
        with open(os.path.join(directory, name), "w") as handle:
            handle.write("# synthetic domtblout\n")
            for line in lines:
                handle.write(line + "\n")
    return directory


def _load(files, cutoff=None):
    folder = _write_folder(files)
    return compleasm.load_hmmsearch_output(folder, CUTOFF if cutoff is None else cutoff)


class TestHmmCandidateRetention(unittest.TestCase):
    def test_two_passing_candidates_are_both_retained(self):
        target_a, protein_a, _ = _target("AAA", "ctgA", 0, 300)
        target_b, protein_b, _ = _target("BBB", "ctgB", 500, 900)
        lines_ab = [
            _domtbl_line(target_a, 80.0, 1, 100),
            _domtbl_line(target_b, 70.0, 1, 90),
        ]
        lines_ba = list(reversed(lines_ab))

        retained_ab = set(_load([("10at2.out", lines_ab)])[0])
        retained_ba = set(_load([("10at2.out", lines_ba)])[0])

        self.assertEqual(retained_ab, {target_a, target_b})
        self.assertEqual(retained_ba, {target_a, target_b})
        self.assertEqual(retained_ab, retained_ba)

    def test_row_order_permutations_are_identical(self):
        specs = [
            ("AAA", "ctgA", 0, 300, 90.0, 1, 80),
            ("BBB", "ctgB", 10, 400, 80.0, 5, 60),
            ("CCC", "ctgC", 20, 500, 70.0, 2, 40),
        ]
        lines = []
        expected = set()
        for reference, contig, start, stop, score, hmm_from, hmm_to in specs:
            target, _, _ = _target(reference, contig, start, stop)
            lines.append(_domtbl_line(target, score, hmm_from, hmm_to))
            expected.add(target)

        permutations = [
            lines,
            list(reversed(lines)),
            [lines[1], lines[2], lines[0]],
            [lines[2], lines[0], lines[1]],
            [lines[2], lines[1], lines[0]],
            [lines[0], lines[2], lines[1]],
        ]
        normalized = []
        for permutation in permutations:
            mappings, qlen, length, tlen, evidence = _load([("10at2.out", permutation)])
            self.assertEqual(set(mappings), expected)
            signature = (
                tuple(mappings),
                tuple(sorted(qlen.items())),
                tuple(sorted(length.items())),
                tuple(sorted(tlen.items())),
                tuple(
                    (
                        item["target_name"],
                        item["hmm_score"],
                        item["qlen"],
                        item["tlen"],
                        item["matched_length"],
                        tuple(item["domains"]),
                    )
                    for item in evidence
                ),
            )
            normalized.append(signature)
        self.assertEqual(len(set(normalized)), 1)

    def test_failing_candidate_is_dropped_even_when_it_appears_first(self):
        target_a, _, _ = _target("AAA", "ctgA", 0, 300)
        target_b, _, _ = _target("BBB", "ctgB", 500, 900)
        orders = [
            [
                _domtbl_line(target_a, 80.0, 1, 50),
                _domtbl_line(target_b, 10.0, 1, 50),
            ],
            [
                _domtbl_line(target_b, 10.0, 1, 50),
                _domtbl_line(target_a, 80.0, 1, 50),
            ],
        ]
        for lines in orders:
            mappings, _, length, _, evidence = _load([("10at2.out", lines)])
            self.assertEqual(set(mappings), {target_a})
            self.assertEqual(set(length), {target_a})
            self.assertEqual([item["target_name"] for item in evidence], [target_a])

    def test_multiple_domains_aggregate_onto_one_candidate(self):
        target, _, _ = _target("AAA", "ctgA", 0, 300)
        # Disjoint domains: (10-1) + (40-20) = 29. A second candidate must not appear.
        lines = [
            _domtbl_line(target, 88.0, 20, 40),
            _domtbl_line(target, 88.0, 1, 10),
        ]
        mappings, _, length, tlen, evidence = _load([("10at2.out", lines)])
        self.assertEqual(mappings, [target])
        self.assertEqual(length[target], 29)
        self.assertEqual(tlen[target], 300)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0]["domains"], [(1, 10), (20, 40)])

        overlap = [
            _domtbl_line(target, 88.0, 1, 30),
            _domtbl_line(target, 88.0, 20, 40),
        ]
        _, _, overlap_length, _, overlap_evidence = _load([("10at2.out", overlap)])
        # Union span added as to-from, matching the upstream interval arithmetic: 29 + (40-30) = 39.
        self.assertEqual(overlap_length[target], 39)
        self.assertEqual(len(overlap_evidence), 1)
        self.assertEqual(overlap_evidence[0]["domains"], [(1, 30), (20, 40)])

    def test_same_genomic_locus_keeps_both_reference_proteins(self):
        target_a, protein_a, location = _target("AAA", "ctg1", 15, 400)
        target_b, protein_b, location_b = _target("BBB", "ctg1", 15, 400)
        self.assertEqual(location, location_b)
        lines = [
            _domtbl_line(target_a, 91.0, 1, 40),
            _domtbl_line(target_b, 77.0, 3, 33),
        ]
        mappings, _, length, _, evidence = _load([("10at2.out", list(reversed(lines)))])
        self.assertEqual(set(mappings), {target_a, target_b})
        self.assertEqual(set(length), {target_a, target_b})
        by_protein = {item["protein_name"]: item for item in evidence}
        self.assertEqual(set(by_protein), {protein_a, protein_b})
        self.assertEqual(by_protein[protein_a]["location"], location)
        self.assertEqual(by_protein[protein_b]["location"], location)
        self.assertNotEqual(by_protein[protein_a]["matched_length"], by_protein[protein_b]["matched_length"])

    def test_distinct_loci_stay_distinct(self):
        target_a, _, location_a = _target("AAA", "ctg1", 0, 100)
        target_b, _, location_b = _target("BBB", "ctg2", 1000, 1500)
        mappings, _, length, _, evidence = _load([
            ("10at2.out", [
                _domtbl_line(target_b, 60.0, 1, 20),
                _domtbl_line(target_a, 60.0, 1, 25),
            ])
        ])
        self.assertEqual(set(mappings), {target_a, target_b})
        self.assertNotEqual(location_a, location_b)
        locations = {item["location"] for item in evidence}
        self.assertEqual(locations, {location_a, location_b})
        self.assertEqual(length[target_a], 24)
        self.assertEqual(length[target_b], 19)

    def test_three_passing_candidates_survive_every_order(self):
        targets = []
        lines = []
        for index, reference in enumerate(("AAA", "BBB", "CCC")):
            target, _, _ = _target(reference, "ctg{}".format(index), index * 100, index * 100 + 50)
            targets.append(target)
            lines.append(_domtbl_line(target, 50.0 + index, 1, 10 + index))
        for permutation in (lines, list(reversed(lines)), [lines[2], lines[0], lines[1]]):
            mappings, _, _, _, evidence = _load([("10at2.out", permutation)])
            self.assertEqual(set(mappings), set(targets))
            self.assertEqual(len(evidence), 3)

    def test_directory_order_does_not_change_candidates(self):
        target_a, _, _ = _target("AAA", "ctgA", 0, 80)
        target_b, _, _ = _target("BBB", "ctgB", 200, 280)
        # Same BUSCO split across two files. Both must be kept, and a reversed
        # directory listing must not overwrite or drop either key.
        folder = _write_folder([
            ("a.out", [_domtbl_line(target_a, 55.0, 1, 15)]),
            ("b.out", [_domtbl_line(target_b, 65.0, 2, 18)]),
        ])

        def _run(names):
            with mock.patch.object(compleasm.os, "listdir", return_value=list(names)):
                return compleasm.load_hmmsearch_output(folder, CUTOFF)

        forward = _run(["a.out", "b.out"])
        reverse = _run(["b.out", "a.out"])
        self.assertEqual(set(forward[0]), {target_a, target_b})
        self.assertEqual(forward[0], reverse[0])
        self.assertEqual(forward[2], reverse[2])
        self.assertEqual(
            [item["target_name"] for item in forward[4]],
            [item["target_name"] for item in reverse[4]],
        )


def _miniprot_frame(rows):
    return pd.DataFrame(rows, columns=[
        "Protein_name", "Contig_name", "Contig_Start", "Contig_Stop",
        "Strand", "Rank", "Identity", "Positive", "Score",
    ])


class TestHmmCandidateEvidenceWriter(unittest.TestCase):
    def test_same_locus_keeps_two_reference_rows_with_miniprot_fields(self):
        target_a, protein_a, location = _target("AAA", "ctg1", 15, 400)
        target_b, protein_b, location_b = _target("BBB", "ctg1", 15, 400)
        self.assertEqual(location, location_b)
        _, _, _, _, evidence = _load([
            ("10at2.out", [
                _domtbl_line(target_b, 77.0, 3, 33),
                _domtbl_line(target_a, 91.0, 1, 40),
            ])
        ])
        frame = _miniprot_frame([
            [protein_b, "ctg1", 15, 400, "-", 2, 0.81, 0.9, 40],
            [protein_a, "ctg1", 15, 400, "+", 1, 0.95, 0.97, 80],
        ])
        with tempfile.TemporaryDirectory(prefix="compleasm-evidence-") as directory:
            path = os.path.join(directory, "hmm_candidate_evidence.tsv")
            rows = compleasm.write_hmm_candidate_evidence(path, evidence, frame)
            with open(path) as handle:
                text = handle.read()
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["protein_name"] for row in rows}, {protein_a, protein_b})
        self.assertEqual({(row["contig"], row["genomic_start"], row["genomic_end"]) for row in rows}, {("ctg1", "15", "400")})
        by_protein = {row["protein_name"]: row for row in rows}
        self.assertEqual(by_protein[protein_a]["strand"], "+")
        self.assertEqual(by_protein[protein_a]["miniprot_rank"], "1")
        self.assertEqual(by_protein[protein_a]["miniprot_identity"], "0.950000")
        self.assertEqual(by_protein[protein_a]["miniprot_positive"], "0.970000")
        self.assertEqual(by_protein[protein_a]["miniprot_score"], "80")
        self.assertEqual(by_protein[protein_a]["hmm_score"], "91.000000")
        self.assertEqual(by_protein[protein_a]["locus_id"], location)
        self.assertIn("compleasm_hmm_evidence_schema_version=1", text)
        self.assertTrue(text.strip().splitlines()[-1].startswith("1\t"))

    def test_writer_bytes_do_not_depend_on_input_order(self):
        target_a, protein_a, _ = _target("AAA", "ctgA", 0, 300)
        target_b, protein_b, _ = _target("BBB", "ctgB", 500, 900)
        forward, _, _, _, evidence_ab = _load([
            ("b.out", [_domtbl_line(target_b, 70.0, 1, 90)]),
            ("a.out", [_domtbl_line(target_a, 80.0, 1, 100)]),
        ])
        reverse, _, _, _, evidence_ba = _load([
            ("a.out", [_domtbl_line(target_a, 80.0, 1, 100)]),
            ("b.out", [_domtbl_line(target_b, 70.0, 1, 90)]),
        ])
        self.assertEqual(forward, reverse)
        frame = _miniprot_frame([
            [protein_b, "ctgB", 500, 900, "-", 1, 0.7, 0.8, 12],
            [protein_a, "ctgA", 0, 300, "+", 1, 0.6, 0.7, 11],
        ])
        with tempfile.TemporaryDirectory(prefix="compleasm-evidence-") as directory:
            path_ab = os.path.join(directory, "ab.tsv")
            path_ba = os.path.join(directory, "ba.tsv")
            compleasm.write_hmm_candidate_evidence(path_ab, list(reversed(evidence_ab)), frame.iloc[::-1])
            compleasm.write_hmm_candidate_evidence(path_ba, evidence_ba, frame)
            with open(path_ab, "rb") as handle_ab, open(path_ba, "rb") as handle_ba:
                self.assertEqual(handle_ab.read(), handle_ba.read())


if __name__ == "__main__":
    unittest.main()
