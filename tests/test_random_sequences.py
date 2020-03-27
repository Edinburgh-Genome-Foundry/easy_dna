from easy_dna import random_dna_sequence, random_protein_sequence


def test_random_protein_sequence():
    n = 10
    random_protein = random_protein_sequence(n, seed=2)
    assert random_protein[0] == "M"
    assert random_protein[-1] == "*"
    assert len(random_protein) == n
    aa_set = set("ACEDGFIHKLNQPSRTWVY" + "M*")
    assert len(aa_set | set(random_protein)) == len(aa_set)
