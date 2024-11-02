import os

from easy_dna import (
    annotate_record,
    sequence_to_biopython_record,
    anonymized_record,
    censor_record,
    censor_genbank,
)

RECORDS_FOLDER = os.path.join("tests", "data")


def test_anonymized_record():
    record = sequence_to_biopython_record(
        "ACGTGCGATGGGATTATTTCCAAC", id="test id", name="test name"
    )
    annotate_record(record, location="full", feature_type="test_feature")
    record = anonymized_record(record)
    assert record.seq == "ACGTGCGATGGGATTATTTCCAAC"
    assert record.id == "anonymized"
    assert record.name == "Exported"
    assert record.features[0].qualifiers["label"] == "feature_0"


def test_censor_record():
    sequence = "ACGTGCGAT" + "CGTCTC" + "ATGC" + "A" * 20
    record = sequence_to_biopython_record(sequence, id="test id", name="test name")
    annotate_record(record, location="full", feature_type="test_feature")
    record = censor_record(
        record, keep_topology=True, anonymise_features=True, preserve_sites=["BsmBI"]
    )
    assert record.seq != sequence
    assert record.id == "censored"
    assert record.name == "Exported"
    assert record.features[0].qualifiers["label"] == "feature_0"
    assert "CGTCTC" in record.seq


def test_censor_genbank(tmpdir):
    filename = os.path.join(RECORDS_FOLDER, "test.gb")
    target = os.path.join(tmpdir, "censored.gb")
    censor_genbank(filename, target)
    assert os.path.isfile(target)
