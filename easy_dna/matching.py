import itertools
from .biotables import NUCLEOTIDE_TO_REGEXPR, IUPAC_NOTATION
from .biological_operations import reverse_complement
from .record_operations import sequence_to_atgc_string


def dna_pattern_to_regexpr(dna_pattern):
    """Return a regular expression pattern for the provided DNA pattern.

    For instance ``dna_pattern_to_regexpr('ATTNN')`` returns
    ``"ATT[A|T|G|C][A|T|G|C]"``.
    """
    return "".join([NUCLEOTIDE_TO_REGEXPR[n] for n in dna_pattern])


def all_iupac_variants(iupac_sequence):
    """Return all unambiguous possible versions of the given sequence.

    Examples
    ========

    >>> all_iupac_variants('ATN')
    >>> ['ATA', 'ATC', 'ATG', 'ATT']
    """
    return [
        "".join(nucleotides)
        for nucleotides in itertools.product(
            *[IUPAC_NOTATION[n] for n in iupac_sequence]
        )
    ]


def find_index(seq, pattern):
    """Return the index of the first match of the pattern in seq.

    Returns -1 if there is no match.
    """

    pattern = sequence_to_atgc_string(pattern)
    if hasattr(seq, "seq"):
        # Seq is a SeqRecord
        return seq.seq.find(pattern)
    else:
        return seq.find(pattern)


def find_occurence(seq, pattern, strand="both"):
    if strand == 1:
        position = find_index(seq, pattern)
        if position == -1:
            return None
        else:
            return (position, position + len(pattern), 1)
    elif strand == -1:
        rev = reverse_complement(seq)
        occurrence = find_occurence(rev, pattern, strand=1)
        if occurrence is None:
            return None
        L = len(seq)
        return (L - occurrence[1], L - occurrence[0], -1)
    elif strand == "both":
        occurrence = find_occurence(seq, pattern, strand=1)
        if occurrence is None:
            occurrence = find_occurence(seq, pattern, strand=-1)
        return occurrence
