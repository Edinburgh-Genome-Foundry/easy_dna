from copy import deepcopy
from io import BytesIO, StringIO
import itertools
import os
import re

import pandas

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import crazydoc
import flametree
from snapgene_reader import snapgene_file_to_seqrecord

from .record_operations import sequence_to_biopython_record

strands = ["   ", "ss-", "ds-", "ms-"]
moleculetypes = [
    "NA    ",
    "DNA   ",
    "RNA   ",
    "tRNA  ",
    "rRNA  ",
    "mRNA  ",
    "uRNA  ",
    "cRNA  ",
]

strand_molecules = [
    "".join(entry).strip(" ") for entry in itertools.product(strands, moleculetypes)
]

genbank_reference = {
    # From the Genbank file format definition (release 260).
    # 0-based slicing indices below.
    # name: (start, end, [options])
    "locus": (0, 5, ["LOCUS"]),
    "name": (12, 28, None),
    "length": (29, 40, None),
    "bp": (41, 43, ["bp"]),
    "strandedness": (44, 47, strands),
    "moleculetype": (
        47,
        53,
        moleculetypes,
    ),
    # No topology info (spaces) also allowed, for flexibility:
    "moleculetopology": (55, 63, ["linear  ", "circular", "        "]),
    "divisioncode": (64, 67, None),
    "updatedate": (68, 79, None),
    # Spaces:
    "space1": (5, 12, "       "),
    "space2": (28, 29, " "),
    "space3": (40, 41, " "),
    "space4": (43, 44, " "),
    "space5": (53, 55, "  "),
    "space6": (63, 64, " "),
    "space7": (67, 68, " "),
}

genbank_months = [
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
]

locus_guide = """+---+       +--------------+ +---------+ +- +--+----+  +------+ +-- +---------+
1   5       13            28 30       40 42 45 48  53  56    63 65  69       79"""


# Guidelines:
# Use only alphanumeric and underscore characters in the sequence name.
# Sequence names must be maximum 20 character long.
# Set SeqRecord.id to the sequence name, and .name="Exported" (this is used in exports).
# Genbank format must be: LOCUS = Exported, filename = sequence name.
# 'Exported' as name is within required length and also used by other software.
def load_record(
    filename, record_id="filename", adapt_id=True, upperize=False, id_cutoff=20
):
    """Load a FASTA/Genbank/Snapgene file as a Biopython record.

    Parameters
    ----------
    filename : str
        Path to the sequence file.
    record_id : str
        ID of the record ("filename": use the file name; "original": keep the record's
        original ID (defaults to file name if the record has no ID).
    adapt_id : bool
        If True, convert ID to alphanumeric and underscore.
    upperize : bool
        If True, the record's sequence will converted to uppercase.
    id_cutoff : int, optional
        If the ID is longer than this value, it will get truncated at this cutoff to
        conform to guidelines and Genbank name limit. Use `None` for no cutoff.
    """
    # FILEFORMAT
    if filename.lower().endswith((".gb", ".gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith((".fa", ".fasta")):
        record = SeqIO.read(filename, "fasta")
    elif filename.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(filename)
    else:
        raise ValueError("Unknown format for file: %s" % filename)
    # LETTERCASE
    if upperize:
        record = record.upper()
    # NAME
    record.name = "Exported"  # see guidelines above
    if record_id == "original":
        if record.id in [None, "", "<unknown id>", ".", " "]:
            record.id = os.path.splitext(os.path.basename(filename))[0]
    elif record_id == "filename":
        record.id = os.path.splitext(os.path.basename(filename))[0]
    else:  # to catch if a deprecated specification is used
        ValueError("`record_id` must be one of `filename` or `original`!")
    # LENGTH
    if id_cutoff is not None:
        record.id = record.id[:id_cutoff]
    # CHARACTERS
    if adapt_id:
        record.id = "".join([char if char.isalnum() else "_" for char in record.id])

    return record


def write_record(record, target, id_cutoff=20, adapt_id=True):
    """Write a DNA record as Genbank or FASTA via Biopython, with fixes.

    Parameters
    ----------
    record : SeqRecord
        Biopython SeqRecord.
    target : str or StringIO
        Filepath, string, or StringIO instance. Desired sequence format is inferred from
        the ending. If it's a directory, it uses the ID as filename and exports Genbank.
    id_cutoff : int, optional
        If the ID is longer than this value, it will get truncated at this cutoff to
        conform to guidelines and Genbank name limit. Use `None` for no cutoff.
    adapt_id : bool
        If True, convert ID to alphanumeric and underscore.
    """
    record = deepcopy(record)
    record.name = "Exported"  # This is used as LOCUS
    if id_cutoff is not None:
        record.id = record.id[:id_cutoff]
    if adapt_id:
        record.id = "".join([char if char.isalnum() else "_" for char in record.id])
    if has_dna_alphabet:  # for Biopython <1.78
        if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
            record.seq.alphabet = DNAAlphabet()
    record.annotations["molecule_type"] = "DNA"

    if hasattr(target, "open"):
        target = target.open("w")
    else:
        if type(target) is StringIO:
            fmt = "genbank"
        elif target.lower().endswith((".gb", ".gbk")):
            fmt = "genbank"
        elif target.lower().endswith((".fa", ".fasta")):
            fmt = "fasta"
        else:  # directory
            target = os.path.join(target, record.id + ".gb")
            fmt = "genbank"
    SeqIO.write(record, target, fmt)


def string_to_record(string):
    """Convert a string of a FASTA, Genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    # print("============", len(matches.groups()[0]), len(string))
    # print (matches.groups()[0] == string)
    if (matches is not None) and (matches.groups()[0] == string):
        if has_dna_alphabet:  # Biopython <1.78
            sequence = Seq(string, alphabet=DNAAlphabet())
        else:
            sequence = Seq(string)
        seqrecord = SeqRecord(sequence)
        seqrecord.annotations["molecule_type"] = "DNA"

        return seqrecord, "ATGC"

    for fmt in ("fasta", "genbank"):
        try:
            stringio = StringIO(string)
            records = list(SeqIO.parse(stringio, fmt))
            if len(records) > 0:
                return (records, fmt)
        except Exception:
            pass
    try:
        record = snapgene_file_to_seqrecord(filecontent=StringIO(string))
        return record
    except Exception:
        pass
    raise ValueError("Invalid sequence format")


def spreadsheet_file_to_dataframe(filepath, header="infer"):
    """Load a CSV or EXCEL spreadsheet as a Pandas dataframe."""
    name = filepath._name if hasattr(filepath, "_name") else filepath
    if name.endswith(".csv"):
        return pandas.read_csv(filepath, header=header)
    else:
        return pandas.read_excel(filepath, header=header)


def records_from_zip_file(zip_file):
    """Return all fasta/genbank/snapgene in a zip as Biopython records."""
    zip_file = flametree.file_tree(zip_file)
    records = []
    for f in zip_file._all_files:
        ext = f._extension.lower()
        if ext in ["gb", "gbk", "fa", "dna"]:
            try:
                new_records, fmt = string_to_record(f.read())
            except Exception:
                content_stream = BytesIO(f.read("rb"))
                try:
                    record = snapgene_file_to_seqrecord(fileobject=content_stream)
                    new_records, _ = [record], "snapgene"
                except Exception:
                    try:
                        parser = crazydoc.CrazydocParser(
                            ["highlight_color", "bold", "underline"]
                        )
                        new_records = parser.parse_doc_file(content_stream)
                        # fmt = "doc"
                    except Exception:
                        raise ValueError("Format not recognized for file " + f._path)

            single_record = len(new_records) == 1
            for i, record in enumerate(new_records):
                name = record.id
                if name in [
                    None,
                    "",
                    "<unknown id>",
                    ".",
                    " ",
                    "<unknown name>",
                ]:
                    number = "" if single_record else ("%04d" % i)
                    name = f._name_no_extension.replace(" ", "_") + number
                record.id = name
                record.name = name
                record.file_name = f._name_no_extension
            records += new_records
    return records


def records_from_file(filepath):
    """Autodetect file format and load Biopython records from it."""

    with open(filepath, "rb") as f:
        content = f.read()
    try:
        records, fmt = string_to_record(content.decode("utf-8"))
    except Exception:
        try:
            record = snapgene_file_to_seqrecord(fileobject=BytesIO(content))
            records, fmt = [record], "snapgene"
        except Exception:
            try:
                parser = crazydoc.CrazydocParser(
                    ["highlight_color", "bold", "underline"]
                )
                records = parser.parse_doc_file(BytesIO(content))
                fmt = "doc"
            except Exception:
                try:
                    df = spreadsheet_file_to_dataframe(filepath, header=None)
                    records = [
                        sequence_to_biopython_record(sequence=seq, id=name, name=name)
                        for name, seq in df.values
                    ]
                    fmt = "spreadsheet"
                except Exception:
                    raise ValueError("Format not recognized for file " + filepath)
    if not isinstance(records, list):
        records = [records]
    return records, fmt


def record_to_formated_string(record, remove_descr=False):
    """Return a string with the content of a Genbank file."""
    if remove_descr:
        record = deepcopy(record)
        if isinstance(record, (list, tuple)):
            for r in record:
                r.description = ""
        else:
            record.description = ""
    fileobject = StringIO()
    write_record(record=record, target=fileobject)

    return fileobject.getvalue().encode("utf-8")


def records_from_data_files(filepaths=None, folder=None):
    """Automatically convert files or a folder's content to Biopython records."""
    if folder is not None:
        filepaths = [f._path for f in flametree.file_tree(folder)._all_files]
    records = []
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        if filename.lower().endswith("zip"):
            records += records_from_zip_file(filepath)
            continue
        recs, fmt = records_from_file(filepath)
        single_record = len(recs) == 1
        for i, record in enumerate(recs):
            name_no_extension = "".join(filename.split(".")[:-1])
            name = name_no_extension + ("" if single_record else ("%04d" % i))
            name = name.replace(" ", "_")
            UNKNOWN_IDS = [
                "None",
                "",
                "<unknown id>",
                ".",
                "EXPORTED",
                "<unknown name>",
                "Exported",
            ]
            if has_dna_alphabet:  # Biopython <1.78
                record.seq.alphabet = DNAAlphabet()
            record.annotations["molecule_type"] = "DNA"

            # Sorry for this parts, it took a lot of "whatever works".
            # keep your part names under 20c and pointless, and everything
            # will be good
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            if str(record.name).strip() in UNKNOWN_IDS:
                record.name = name
            record.file_name = name_no_extension
        records += recs
    return records


def is_genbank_standard(filepath):
    """Check the LOCUS line of a Genbank file."""
    std_locus_line_len = 79  # as per Genbank definition
    is_correct = True
    expected_entries = 8  # expected, as per Genbank definition
    # Strandedness and Molecule Type is handled as one entry.

    with open(filepath) as genbank:
        first_line = genbank.readline().strip("\n")

    entries = re.split(" +", first_line)
    if len(entries) == expected_entries:
        topology = entries[-3]
    # Topology can be empty (spaces), in that case we expect one fewer column:
    elif len(entries) == expected_entries - 1:
        topology = "        "
    else:  # bad first line, quit
        print("Error: malformatted locus line: missing or too many entries")
        print(first_line, locus_guide, sep="\n")
        is_correct = False

        return is_correct

    entry_dict = {
        "locus": entries[0],
        "name": entries[1],
        "length": entries[2],
        "bp": entries[3],
        "strand_molecules": entries[4],
        "moleculetopology": topology,
        "divisioncode": entries[-2],
        "updatedate": entries[-1],
    }

    if len(first_line) == std_locus_line_len:
        # Look up each entry:
        for key, value in genbank_reference.items():
            if value[2] is not None:
                if first_line[value[0] : value[1]] not in value[2]:
                    print("Error:", key, "=", first_line[value[0] : value[1]])
                    print("       It should be one of ", value[2])
                    is_correct = False
    else:  # incorrect length, but we try and parse using whitespace
        print(
            "Error: LOCUS line is %d characters long! (standard length = %d)"
            % (
                len(first_line),
                std_locus_line_len,
            )
        )
        is_correct = False

        # This checks the same as the block above, except on whitespace-split entries:
        for entry, value in entry_dict.items():
            if entry in genbank_reference.keys():
                if genbank_reference[entry][2] is not None:
                    if value not in genbank_reference[entry][2]:
                        print("Error:", entry, "=", value)
                        print(
                            "       It should be one of ", genbank_reference[entry][2]
                        )
        # We split by whitespace, therefore this entry won't be in the Genbank ref
        if entry_dict["strand_molecules"] not in strand_molecules:
            print(
                "Error: incorrect strandness+molecule type = ",
                entry_dict["strand_molecules"],
            )

    # Inspecting entries in detail:
    if not entry_dict["length"].isnumeric():
        is_correct = False
        print("Error: non-numeric characters in length:", entries[2])
    # Date should be dd-MMM-yyyy (15-NOV-2024)
    date_values = entry_dict["updatedate"].split("-")
    if not date_values[0].isnumeric():  # ignores non-valid day number (e.g. 00)
        is_correct = False
        print("Error: non-numeric characters in day:", date_values[0])
    if date_values[1] not in genbank_months:
        is_correct = False
        print("Error: invalid month:", date_values[1])
    if not date_values[2].isnumeric():  # ignores non-valid year number
        is_correct = False
        print("Error: non-numeric characters in year:", date_values[2])

    if not is_correct:
        print(first_line, locus_guide, sep="\n")

    return is_correct
