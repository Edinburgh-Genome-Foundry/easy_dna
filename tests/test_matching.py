from easy_dna import dna_pattern_to_regexpr, all_iupac_variants
from easy_dna.matching import find_index


def test_dna_pattern_to_regexpr():
    assert dna_pattern_to_regexpr("ATW") == "AT[ATW]"


def test_all_iupac_variants():
    assert sorted(all_iupac_variants("ATN")) == ["ATA", "ATC", "ATG", "ATT"]


def test_find_index():
    assert find_index("ATCGATTCGATC", "CGA") == 2
