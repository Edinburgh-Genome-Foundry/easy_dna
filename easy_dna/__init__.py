from .biological_operations import (
    complement,
    reverse_complement,
    translate,
    reverse_translate,
)

from .biotables import (
    CODONS_SEQUENCES,
    CODONS_SYNONYMS,
    CODONS_TRANSLATIONS,
    COMPLEMENTS,
)

from .io import load_record, write_record, records_from_data_files

from .matching import dna_pattern_to_regexpr, all_iupac_variants

from .random_sequences import random_dna_sequence, random_protein_sequence

from .enzymes import list_common_enzymes

from .record_operations import (
    annotate_record,
    sequence_to_biopython_record,
    record_with_different_sequence,
    anonymized_record,
    censor_record,
    censor_genbank,
)

from .modifications import (
    replace_segment,
    replace_occurence,
    insert_segment,
    delete_segment,
    delete_nucleotides,
    reverse_segment,
    cut_and_paste_segment,
    copy_and_paste_segment,
    swap_segments,
)

from .extractor import extract_from_input

from .version import __version__

__all__ = [
    "complement",
    "reverse_complement",
    "translate",
    "reverse_translate",
    "CODONS_SEQUENCES",
    "CODONS_SYNONYMS",
    "CODONS_TRANSLATIONS",
    "COMPLEMENTS",
    "load_record",
    "write_record",
    "anonymized_record",
    "censor_record",
    "censor_genbank",
    "records_from_data_files",
    "dna_pattern_to_regexpr",
    "all_iupac_variants",
    "random_dna_sequence",
    "random_protein_sequence",
    "list_common_enzymes",
    "annotate_record",
    "sequence_to_biopython_record",
    "record_with_different_sequence",
    "replace_segment",
    "replace_occurence",
    "insert_segment",
    "delete_segment",
    "delete_nucleotides",
    "reverse_segment",
    "cut_and_paste_segment",
    "copy_and_paste_segment",
    "swap_segments",
    "extract_from_input",
]
