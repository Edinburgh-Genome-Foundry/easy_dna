from easy_dna import (
    replace_segment,
    replace_occurence,
    insert_segment,
    delete_segment,
    delete_nucleotides,
    reverse_segment,
    cut_and_paste_segment,
    copy_and_paste_segment,
    swap_segments,
)


def test_replace_segment():
    assert replace_segment("ATGCATGCTACGGGCT", 2, 5, "AAAAAA") == "ATAAAAAATGCTACGGGCT"


def test_insert_segment():
    assert insert_segment("ATGCGGGCT", 2, "AAAAAA") == "ATAAAAAAGCGGGCT"


def test_delete_segment():
    assert delete_segment("ATGCATGCTACGGGCT", 2, 5) == "ATTGCTACGGGCT"


def test_delete_nucleotides():
    assert delete_nucleotides("AAAAAAGCGGGCT", 0, 6) == "GCGGGCT"


def test_reverse_segment():
    assert reverse_segment("AAAAGTTCCAAAA", 4, 9) == "AAAAGGAACAAAA"


def test_cut_and_paste_segment():
    assert cut_and_paste_segment("AAAAGTTCCAAAA", 4, 9, 0) == "GTTCCAAAAAAAA"
    assert cut_and_paste_segment("AAAAGTTCCAAAA", 4, 9, 13) == "AAAAAAAAGTTCC"


def test_replace_occurence():
    assert replace_occurence("AAAGTTCCA", "GGAA", "CCCC", strand="both") == "AAAGGGGGA"


def test_swap_segments():
    assert swap_segments("AAAGTTCCA", pos1=(0, 3), pos2=(8, 9)) == "AGTTCCAAA"


def test_copy_and_paste_segment():
    assert copy_and_paste_segment("AAAGTTCCA", 0, 3, 9) == "AAAGTTCCAAAA"
