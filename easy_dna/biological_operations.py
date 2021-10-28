import numpy as np

from Bio.Seq import Seq

from .biotables import COMPLEMENTS, CODONS_SEQUENCES


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses Biopython for speed.
    """
    if hasattr(dna_sequence, "complement"):
        return dna_sequence.complement()
    if len(dna_sequence) <= 30:
        return "".join([COMPLEMENTS[nuc] for nuc in dna_sequence])
    # This alternative has overhead but is really fast on long sequences
    return str(Seq(dna_sequence).complement())


def reverse_complement(dna_sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``reverse_complement("ATGCCG")`` returns ``"CGGCAT"``.

    Uses Biopython for speed.
    """
    if hasattr(dna_sequence, "reverse_complement"):
        return dna_sequence.reverse_complement()
    return complement(dna_sequence)[::-1]


def reverse_translate(protein_sequence, randomize_codons=False):
    """Return a DNA sequence which translates to the provided protein sequence.

    Note: at the moment, the first valid codon found is used for each
    amino-acid (so it is deterministic but no codon-optimization is done).
    """
    if randomize_codons:
        random_indices = np.random.randint(0, 1000, len(protein_sequence))
        return "".join(
            [
                CODONS_SEQUENCES[aa][random_index % len(CODONS_SEQUENCES[aa])]
                for aa, random_index in zip(protein_sequence, random_indices)
            ]
        )
    return "".join([CODONS_SEQUENCES[aa][0] for aa in protein_sequence])


def translate(dna_sequence, translation_table="Bacterial"):
    """Translate the DNA sequence into an amino-acid sequence "MLKYQT...".

    If ``translation_table`` is the name or number of a NCBI genetic table,
    Biopython will be used. See here for options:

    http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec26

    ``translation_table`` can also be a dictionary of the form
    ``{"ATT": "M", "CTC": "X", etc.}`` for more exotic translation tables.
    """
    if isinstance(translation_table, dict):
        return "".join(
            [
                translation_table[dna_sequence[i : i + 3]]
                for i in range(0, len(dna_sequence), 3)
            ]
        )
    else:
        return str(Seq(dna_sequence).translate(table=translation_table))
