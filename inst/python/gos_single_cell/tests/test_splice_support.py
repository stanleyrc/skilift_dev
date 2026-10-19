"""Unit tests for gos_sc_splice_support.py (synthetic fixtures only)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_splice_support as S  # noqa: E402

EV = [{"gene": "G1", "type": "Splice", "Variant": "c.1-1G>C", "Genome_Location": "17:1000-1000", "seqnames": "17", "n_cells": 10},
      {"gene": "G2", "type": "Splice", "Variant": "c.5+1G>A", "Genome_Location": "3:5000-5000", "seqnames": "3", "n_cells": 4},
      {"gene": "G3", "type": "Missense", "Genome_Location": "1:1-1"}]
SP = {"findings": [{"gene": "G1", "kind": "patient_specific", "severity": "high", "text": "G1 cryptic acceptor", "junction": {"chromosome": "17", "start": 1044, "end": 3000},
                    "numbers": {"n_cells": 11, "reads": 58, "total": 59}},
                   {"gene": "G4", "kind": "known_variant", "severity": "high", "text": "G4 skip", "junction": {"chromosome": "7", "start": 1, "end": 2}, "cells": ["a", "b"]}],
      "events": [{"gene": "G2", "annotated": False, "cells_using": 1, "chromosome": "3", "start": 5010, "end": 9000},
                 {"gene": "G2", "annotated": True, "cells_using": 5, "chromosome": "3", "start": 5010, "end": 9000}]}


class T(unittest.TestCase):
    def test_support_classes(self):
        out = S.patient_support(EV, SP)
        self.assertEqual(out["G1"]["support"], "DNA+RNA")
        self.assertEqual(out["G1"]["rna"][0]["reads"], 58)
        # G2: only a one-cell unannotated event and an annotated one: not confirmed
        self.assertEqual(out["G2"]["support"], "DNA")
        self.assertEqual(out["G4"]["support"], "RNA")
        self.assertEqual(out["G4"]["rna"][0]["n_cells"], 2)
        self.assertNotIn("G3", out)

    def test_far_junction_does_not_confirm(self):
        sp = {"findings": [dict(SP["findings"][0], junction={"chromosome": "17", "start": 5000, "end": 9000})], "events": []}
        self.assertEqual(S.patient_support(EV, sp)["G1"]["support"], "DNA")

    def test_no_rna(self):
        out = S.patient_support(EV, None)
        self.assertEqual({g: v["support"] for g, v in out.items()}, {"G1": "DNA", "G2": "DNA"})


if __name__ == "__main__":
    unittest.main()
