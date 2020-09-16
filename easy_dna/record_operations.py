from copy import deepcopy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from Bio.SeqFeature import SeqFeature, FeatureLocation


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    if hasattr(sequence, "seq"):
        return sequence
    if has_dna_alphabet:  # Biopython <1.78
        sequence = Seq(sequence, alphabet=DNAAlphabet())
    else:
        sequence = Seq(sequence)

    seqrecord = SeqRecord(sequence, id=id, name=name, features=list(features),)
    seqrecord.annotations["molecule_type"] = "DNA"

    return seqrecord


def sequence_to_atgc_string(sequence):
    if isinstance(sequence, str):
        return sequence
    else:
        return str(sequence.seq)


def record_with_different_sequence(record, new_seq):
    """Return a version of the record with the sequence set to new_seq."""
    new_record = deepcopy(record)
    if has_dna_alphabet:  # Biopython <1.78
        sequence = Seq(new_seq, alphabet=DNAAlphabet())
    else:
        sequence = Seq(new_seq)

    new_record.seq = sequence

    return new_record


def annotate_record(
    seqrecord, location="full", feature_type="misc_feature", margin=0, **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The Biopython SeqRecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1).

    feature_type
      The type associated with the feature.

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionary that will be the Biopython feature's `qualifiers` attribute.
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


def anonymized_record(record, record_id="anonymized", label_generator="feature_%d"):
    """Return a record with removed annotations/keywords/features/etc.

    Warning: this does not change the record sequence!

    Parameters
    ----------
    record
      The record to be anonymized.

    record_id
      ID of the new record.

    label_generator
      Recipe to change feature labels. Either "feature_%d" or None (no label)
      of a function (i, feature)=>label.
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
