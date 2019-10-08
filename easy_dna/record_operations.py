from copy import deepcopy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    if hasattr(seq, 'seq'):
        return seq
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=name,
        features=list(features),
    )

def sequence_to_atgc_string(sequence):
    if isinstance(sequence, str):
        return sequence
    else:
        return str(sequence.seq)


def record_with_different_sequence(record, new_seq):
    """Return a version of the record with the sequence set to new_seq"""
    new_record = deepcopy(record)
    new_record.seq = Seq(new_seq, alphabet=DNAAlphabet())
    return new_record


def annotate_record(
    seqrecord,
    location="full",
    feature_type="misc_feature",
    margin=0,
    **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )
