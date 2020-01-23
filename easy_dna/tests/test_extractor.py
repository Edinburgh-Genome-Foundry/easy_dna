import numpy as np
import easy_dna as dna


def test_extract_from_input():
    parts = []
    for i in range(10):
        part_id = "part_%s" % ("ABCDEFGHAB"[i])  # id is nonunique on purpose
        alias = "part_%d" % i  # alias is unique
        part_length = np.random.randint(1000, 1500)
        sequence = dna.random_dna_sequence(part_length)
        record = dna.sequence_to_biopython_record(sequence, id=part_id)
        record.name = part_id
        dna.annotate_record(record, label=part_id, alias=alias)
        parts.append(record)

    constructs = []
    for position_of_last_part in [8, 10]:
        # 8: parts A-H; 10: parts A-H and AB again
        construct_record = sum(parts[1:position_of_last_part], parts[0])
        construct_record.id = "construct_%02d" % (position_of_last_part)
        construct_record.name = construct_record.id
        constructs.append(construct_record)

    processed_report = dna.extract_from_input(construct_list=constructs)
    assert sum(processed_report["has_copy"]) == 8
