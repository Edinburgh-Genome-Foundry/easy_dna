from copy import deepcopy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    if hasattr(sequence, "seq"):
        return sequence
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
            FeatureLocation(int(location[0]), int(location[1]), strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )


def anonymized_record(
    record, record_id="anonymized", label_generator="feature_%d"
):
    """Return a record with removed annotations/keywords/features/etc.

    Warning: this does not change the record sequence!

    Parameters
    ----------
    record
      The record to be anonymized

    record_id
      Id of the new record

    label_generator
      Recipe to change feature labels. Either "feature_%d" or None (no label)
      of a function (i, feature)=>label
    """
    new_record = deepcopy(record)
    new_record.annotations = {
        "molecule_type": "ds-DNA",
        "data_file_division": "   ",
        "keywords": [],
    }
    new_record.id = new_record.name = record_id
    for i, feature in enumerate(new_record.features):
        label = None
        if hasattr(label_generator, "__call__"):
            label = label_generator(i, feature)
        elif isinstance(label_generator, str):
            label = label_generator % i
        feature.qualifiers = {}
        if label is not None:
            feature.qualifiers["label"] = label
    return new_record
