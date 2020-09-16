from Bio.Seq import Seq
from easy_dna import (
    complement,
    reverse_complement,
    translate,
    reverse_translate,
    CODONS_TRANSLATIONS,
)


def test_complement():
    assert set(complement(Seq("A" * 31))) == set("T")
    assert set(complement("A" * 31)) == set("T")


def test_reverse_complement():
    assert reverse_complement(Seq("ATCG")) == str("CGAT")


def test_translate():
    assert translate("AAAAAC") == "KN"
    assert translate("AAAAAC", translation_table=CODONS_TRANSLATIONS) == "KN"


def test_reverse_translate():
    assert reverse_translate("KN") == "AAAAAC"
