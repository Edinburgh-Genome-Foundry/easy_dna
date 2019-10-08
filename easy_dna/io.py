import os
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from snapgene_reader import snapgene_file_to_seqrecord
from copy import deepcopy


def load_record(filename, linear=True, record_id="auto", upperize=False):
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
    record.linear = linear
    if record_id == "auto":
        record_id = record.id
        if record_id in [None, "", "<unknown id>", ".", " "]:
            record_id = os.path.splitext(os.path.basename(filename))[0]
            record.id = record_id
            record.name = record_id.replace(" ", "_")[:20]
        record.id = record_id
    elif record_id is not None:
        record.id = record_id
        record.name = record_id.replace(" ", "_")[:20]
    return record


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.name = record.name[:20]
    if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)
