from easy_dna import dna_pattern_to_regexpr, all_iupac_variants
from easy_dna.matching import find_index, find_occurence


def test_dna_pattern_to_regexpr():
    assert dna_pattern_to_regexpr("ATW") == "AT[ATW]"


def test_all_iupac_variants():
    assert sorted(all_iupac_variants("ATN")) == ["ATA", "ATC", "ATG", "ATT"]


def test_find_index():
    assert find_index("ATCGATTCGATC", "CGA") == 2


def test_find_occurence():
    assert find_occurence('ATCGATTCGATC', 'CGA', strand="both") == (2, 5, 1)
    assert find_occurence('ATCGATTCGATC', 'AAAAA', strand="both") == None
    assert find_occurence('ATCGATTCGATC', 'AATC', strand="both") == (3, 7, -1)
