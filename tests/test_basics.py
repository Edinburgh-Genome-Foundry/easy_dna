from easy_dna import (
    dna_pattern_to_regexpr,
    record_with_different_sequence,
    sequence_to_biopython_record,
    annotate_record,
    list_common_enzymes,
)


def test_record_with_different_sequence():
    record = sequence_to_biopython_record("ATGCATGCATGC")
    annotate_record(record, (0, 5), label="my_label")
    new_record = record_with_different_sequence(record, "GGCCGGCCGGCCGGCC")
    assert len(new_record.features) == 1
    assert new_record.features[0].location == record.features[0].location


def test_list_common_enzymes():
    assert len(list_common_enzymes(min_suppliers=3)) >= 61
    # biopython v1.80: 61 enzymes
    # < v1.79: 63 enzymes
