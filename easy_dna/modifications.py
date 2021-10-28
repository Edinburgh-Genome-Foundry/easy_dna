from .biological_operations import reverse_complement
from .matching import find_index


def replace_segment(seq, start, end, replacement):
    """Return the sequence with ``seq[start:end]`` replaced by ``replacement``."""
    return seq[:start] + replacement + seq[end:]


def insert_segment(seq, pos, inserted):
    """Return the sequence with ``inserted`` inserted, starting at index ``pos``."""
    return seq[:pos] + inserted + seq[pos:]


def delete_segment(seq, start, end):
    """Return the sequence with deleted segment from ``start`` to ``end``."""
    return seq[:start] + seq[end:]


def delete_nucleotides(seq, start, n):
    """Return the sequence with ``n`` deletions from position ``start``."""
    return seq[:start] + seq[start + n :]


def reverse_segment(seq, start, end):
    """Return the sequence with segment ``seq[start:end]`` reverse-complemented."""
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]


def cut_and_paste_segment(seq, start, end, new_start):
    """Move a subsequence by "diff" nucleotides the left or the right."""
    sub = seq[start:end]
    diff = new_start - start
    if diff > 0:
        return seq[:start] + seq[end : end + diff] + sub + seq[end + diff :]
    else:
        return seq[: start + diff] + sub + seq[start + diff : start] + seq[end:]


def replace_occurence(seq, pattern, replacement, strand="both"):
    if strand == -1:
        rev_seq = replace_occurence(
            reverse_complement(seq), pattern, replacement, strand=1
        )
        return reverse_complement(rev_seq)
    elif strand == "both":
        new_seq = replace_occurence(seq, pattern, replacement, strand=1)
        if new_seq is seq:
            return replace_occurence(new_seq, pattern, replacement, strand=-1)
        else:
            return new_seq
    position = find_index(seq, pattern)
    if position == -1:
        return seq
    return replace_segment(seq, position, position + len(pattern), replacement)


def swap_segments(seq, pos1, pos2):
    """Return a new sequence with segments at position ``pos1`` and ``pos2`` swapped.

    ``pos1``, ``pos2`` are both of the form (start1, end1), (start2, end2).
    """
    (start1, end1), (start2, end2) = sorted([pos1, pos2])
    return (
        seq[:start1]
        + seq[start2:end2]
        + seq[end1:start2]
        + seq[start1:end1]
        + seq[end2:]
    )


def copy_and_paste_segment(seq, start, end, new_start):
    """Return the sequence with segment ``[start, end]`` also copied elsewhere,
    starting in ``new_start`."""
    return insert_segment(seq, new_start, seq[start:end])
