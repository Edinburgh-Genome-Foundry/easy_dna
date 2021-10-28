from copy import deepcopy

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Restriction

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False

from .random_sequences import random_dna_sequence
import easy_dna


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
      Recipe to change feature labels. Either ``"feature_%d"`` or ``None`` (no label)
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


def censor_record(
    record,
    record_id="censored",
    label_generator="feature_%d",
    keep_topology=False,
    anonymise_features=True,
    preserve_sites=None,
):
    """Return a record with random sequence and censored annotations/features.

    Useful for creating example files or anonymising sequences for bug reports.


    Parameters
    ----------

    record
      The record to be anonymized.

    record_id
      ID of the new record.

    label_generator
      Recipe to change feature labels. Either ``"feature_%d"`` or ``None`` (no label)
      of a function (i, feature)=>label.

    keep_topology
      Whether to keep the record topology or not.

    anonymise_features
      Whether to replace feature labels and ID/name, or not.

    preserve_sites
      List of enzyme sites to keep. Example: ``["BsmBI", "BsaI"]``. Preserves the
      sequence around cut sites of the specified enzymes.
    """
    # Anonymise
    if anonymise_features:
        new_record = anonymized_record(
            record, record_id=record_id, label_generator=label_generator
        )
    else:
        new_record = deepcopy(record)

    if keep_topology:
        try:
            new_record.annotations["topology"] = record.annotations["topology"]
        except KeyError:  # input may not have topology set
            pass

    # Randomise
    new_seq = random_dna_sequence(
        len(new_record), gc_share=None, probas=None, seed=None
    )

    if preserve_sites:
        restriction_batch = Restriction.RestrictionBatch(preserve_sites)
        # Destroy random new enzyme sites:
        analysis = Restriction.Analysis(restriction_batch, sequence=Seq(new_seq))
        analysis_results = analysis.full()
        for enzyme, hits in analysis_results.items():
            for hit in hits:
                # 10 bp up- and downstream destroys the site whichever strand it is on:
                if hit - 10 < 0:  # handle edge cases
                    start = 0
                    upstream = "A" * hit
                else:
                    start = hit - 10
                    upstream = "A" * 10
                if hit + 10 > len(new_seq):
                    end = len(new_seq)
                    downstream = "A" * (len(new_seq) - hit)
                else:
                    end = hit + 10
                    downstream = "A" * 10
                replacement = upstream + downstream
                new_seq = easy_dna.replace_segment(new_seq, start, end, replacement)

        # Add original sites:
        analysis = Restriction.Analysis(restriction_batch, sequence=record.seq)
        analysis_results = analysis.full()
        original_seq = str(record.seq)
        for enzyme, hits in analysis_results.items():
            for hit in hits:
                # keep 12 bp surrounding the cut site, to capture enzyme site:
                if hit - 12 < 0:  # handle edge cases
                    start = 0
                else:
                    start = hit - 12
                if hit + 12 > len(new_seq):
                    end = len(new_seq)
                else:
                    end = hit + 12
                original_segment = original_seq[start:end]
                new_seq = easy_dna.replace_segment(
                    new_seq, start, end, original_segment
                )

    censored_record = record_with_different_sequence(new_record, new_seq)

    return censored_record


def censor_genbank(filename, target, **censor_params):
    """Load Genbank file and write censored version.


    Parameters
    ----------

    filename
      Path to the file containing the record.

    target
      Path to output genbank file.

    censor_params
      Optional parameters. See ``censor_record()`` for details.
    """
    record = easy_dna.load_record(filename)
    censored = easy_dna.censor_record(record, **censor_params)
    easy_dna.write_record(censored, target)
