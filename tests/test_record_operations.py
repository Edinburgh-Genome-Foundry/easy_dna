from easy_dna import (
    annotate_record,
    sequence_to_biopython_record,
    record_with_different_sequence,
    anonymized_record,
    censor_record,
)


def test_anonymized_record():
    record = sequence_to_biopython_record(
        "ACGTGCGATGGGATTATTTCCAAC", id="test id", name="test name"
    )
    annotate_record(record, location="full", feature_type="test_feature")
    record = anonymized_record(record)
    assert record.seq == "ACGTGCGATGGGATTATTTCCAAC"
    assert record.id == "anonymized"
    assert record.name == "anonymized"
    assert record.features[0].qualifiers["label"] == "feature_0"


def test_censor_record():
    record = sequence_to_biopython_record(
        "ACGTGCGATGGGATTATTTCCAAC", id="test id", name="test name"
    )
    annotate_record(record, location="full", feature_type="test_feature")
    record = censor_record(record, keep_topology=True)
    assert record.seq != "ACGTGCGATGGGATTATTTCCAAC"
    assert record.id == "censored"
    assert record.name == "censored"
    assert record.features[0].qualifiers["label"] == "feature_0"
