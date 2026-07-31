import gzip
import tempfile
import unittest
from pathlib import Path

from lafa_main import entry_id, load_ontology, propagate_scores, read_fasta, write_predictions


class LafaContractTests(unittest.TestCase):
    def test_fasta_ids_match_lafa_accessions(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "query.fasta"
            path.write_text(">sp|P12345|NAME description\nACDX\n>plain-id\nGG\n")
            self.assertEqual(list(read_fasta(path)), [("P12345", "ACDX"), ("plain-id", "GG")])
        self.assertEqual(entry_id("tr|Q9XYZ1|NAME"), "Q9XYZ1")

    def test_ontology_filter_and_propagation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "go.obo"
            path.write_text(
                "[Term]\nid: GO:0000001\n\n"
                "[Term]\nid: GO:0000002\nis_a: GO:0000001 ! parent\n\n"
                "[Term]\nid: GO:0000003\nis_obsolete: true\n"
            )
            active, parents = load_ontology(path)
            result = propagate_scores(["GO:0000002", "GO:0000003"], [0.7, 0.9], active, parents)
            self.assertEqual(result, {"GO:0000002": 0.7, "GO:0000001": 0.7})

    def test_headerless_gzip_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.tsv.gz"
            write_predictions(path, [("P1", "GO:0000001", 0.125)])
            with gzip.open(path, "rt") as handle:
                self.assertEqual(handle.read(), "P1\tGO:0000001\t0.125000\n")


if __name__ == "__main__":
    unittest.main()
