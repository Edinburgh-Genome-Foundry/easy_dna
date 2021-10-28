from copy import deepcopy
from io import BytesIO, StringIO
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


def load_record(filename, record_id="auto", upperize=False, id_cutoff=20):
    """Load a Fasta/Genbank/Snapgene file as a Biopython record.

    Parameters
    ==========
    filename
      Path to the file containing the record.

    record_id
      Id of the record (leave to "auto" to keep the record's original Id, which
      will default to the file name if the record has no Id).

    upperize
      If true, the record's sequence will be upperized.

    id_cutoff
      If the Id is read from a filename, it will get truncated at this cutoff
      to avoid errors at report write time.
    """
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filename, "fasta")
    elif filename.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(filename)
    else:
        raise ValueError("Unknown format for file: %s" % filename)
    if upperize:
        record = record.upper()
    if record_id == "auto":
        record_id = record.id
        if record_id in [None, "", "<unknown id>", ".", " "]:
            record_id = os.path.splitext(os.path.basename(filename))[0]
            record.id = record_id
            record.name = record_id.replace(" ", "_")[:id_cutoff]
        record.id = record_id
    elif record_id is not None:
        record.id = record_id
        record.name = record_id.replace(" ", "_")[:id_cutoff]
    return record


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes."""
    record = deepcopy(record)
    record.name = record.name[:20]
    if has_dna_alphabet:  # Biopython <1.78
        if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
            record.seq.alphabet = DNAAlphabet()
    record.annotations["molecule_type"] = "DNA"

    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)


def string_to_record(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

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


def record_to_formated_string(record, fmt="genbank", remove_descr=False):
    """Return a string with the content of a FASTA/GENBANK file."""
    if remove_descr:
        record = deepcopy(record)
        if isinstance(record, (list, tuple)):
            for r in record:
                r.description = ""
        else:
            record.description = ""
    fileobject = StringIO()
    write_record(record, fileobject, fmt)
    return fileobject.getvalue().encode("utf-8")


def records_from_data_files(filepaths=None, folder=None):
    """Automatically convert files or a folder's content to Biopython records.
    """
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
