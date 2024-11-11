import os

from pandas import DataFrame

from Bio.SeqRecord import SeqRecord

from easy_dna import load_record, sequence_to_biopython_record
from easy_dna.io import (
    string_to_record,
    spreadsheet_file_to_dataframe,
    records_from_zip_file,
    records_from_file,
    record_to_formated_string,
    records_from_data_files,
    is_genbank_standard,
)

RECORDS_FOLDER = os.path.join("tests", "data")
CSV_FOLDER = os.path.join("easy_dna", "data")


def test_load_record():
    record = load_record(os.path.join(RECORDS_FOLDER, "test.gb"), upperize=True)
    assert isinstance(record, SeqRecord)


def test_string_to_record():
    records = string_to_record(">test\nATGC\n")
    assert len(records[0]) == 1  # list of records
    assert isinstance(records[0][0], SeqRecord)
    assert records[1] == "fasta"


def test_spreadsheet_file_to_dataframe():
    df = spreadsheet_file_to_dataframe(os.path.join(CSV_FOLDER, "complements.csv"))
    assert isinstance(df, DataFrame)


def test_records_from_zip_file():
    records = records_from_zip_file(os.path.join(RECORDS_FOLDER, "test.zip"))
    assert len(records) == 1
    assert isinstance(records[0], SeqRecord)


def test_records_from_file():
    record, fmt = records_from_file(os.path.join(RECORDS_FOLDER, "test.gb"))
    assert isinstance(record[0], SeqRecord)
    assert fmt == "genbank"


def test_record_to_formated_string():
    record = sequence_to_biopython_record(
        "ACGTGCGATGGGATTATTTCCAAC", id="test id", name="test name"
    )

    record_string = record_to_formated_string(record, remove_descr=True)
    assert isinstance(record_string, bytes)


def test_records_from_data_files():
    records = records_from_data_files(folder=RECORDS_FOLDER)
    for record in records:
        assert isinstance(record, SeqRecord)


def test_is_genbank_standard():
    assert is_genbank_standard(os.path.join(RECORDS_FOLDER, "test.gb")) == True
