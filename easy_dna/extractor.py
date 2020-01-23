import pandas as pd
import flametree
from .io import load_record, records_from_data_files, write_record


def extract_from_input(
    file=None, directory=None, construct_list=None, direct_sense=False
):
    """Extract features from input and save them in separate files.

    Parameters
    ==========

    file
      Input sequence file (Genbank).

    directory
      Directory name containing input sequence files.

    construct_list
      A list of SeqRecords.

    direct_sense
      If True: make antisense features into direct-sense in the exported files.
    """
    if construct_list:
        records_dict = dict()
        for input_record in construct_list:
            records = extract_features(input_record)
            # potential clash if names are shared:
            key = input_record.name[0:20]  # GenBank format hard limit for name
            records_dict[key] = records
    else:
        records_dict = run_extraction(
            file=file, directory=directory, direct_sense=False
        )

    for k, v in records_dict.items():
        write_records(k, records_dict)

    parts_report = make_part_dict(records_dict)
    processed_report = process_report(parts_report[1])

    common_parts_dict = parts_report[0]
    common_parts_write = dict()
    common_parts_write["common_parts"] = list(common_parts_dict.values())

    write_records("common_parts", common_parts_write, add_prefix=False)

    return processed_report


def run_extraction(file=None, directory=None, direct_sense=False):
    """Run extract_features() on a Genbank file or directory of files.
    """
    genbank_id_limit = 20  # GenBank format hard limit for name
    if file:
        input_record = load_record(
            file, record_id="auto", upperize=False, id_cutoff=genbank_id_limit
        )
        all_input_records = [input_record]
    elif directory:
        all_input_records = records_from_data_files(filepaths=None, folder=directory)
    else:
        raise TypeError("Specify one of 'file' or 'directory'.")

    records_dict = dict()

    for input_record in all_input_records:
        records = extract_features(input_record)

        # potential clash if names are shared:
        key = input_record.name[0:genbank_id_limit]
        records_dict[key] = records

    return records_dict


def compute_id(seq_feature):
    """Computes an id for extract_features().

    The id can be used as a SeqRecord/Genbank name or id.
    """
    label_fields = [
        "name",
        "gene",
        "label",
        "product",
        "source",
        "note",
    ]

    label = "no_label"
    for key in label_fields:
        if key in seq_feature.qualifiers and len(seq_feature.qualifiers[key]):
            label = seq_feature.qualifiers[key]
            if isinstance(label, list):
                label = "_".join(label)
            break

    genbank_limit = 20  # GenBank format hard limit for name
    label = label[0:genbank_limit]
    label = label.replace(" ", "_").replace("'", "p")
    return label


def extract_features(seq_record, direct_sense=True):
    """Extract all features from a SeqRecord.

    Return the features as a list of SeqRecords.
    """
    all_features_in_direct_sense = direct_sense
    records = []

    # Check whether the feature is fully contained within any other feature:
    for i, seq_feature in enumerate(seq_record.features):
        start = seq_feature.location.start
        end = seq_feature.location.end

        contained_feature = False
        for j, other_seq in enumerate(
            seq_record.features[:i] + seq_record.features[(i + 1) :]
        ):  # exclude current feature
            if start in other_seq.location and end in other_seq.location:
                contained_feature = True
                break

        if contained_feature:
            continue

        if all_features_in_direct_sense:
            new_record = seq_feature.extract(seq_record)
        else:
            new_record = seq_record[
                seq_feature.location.start : seq_feature.location.end
            ]

        new_record.annotations = {"topology": "linear"}
        new_record.description = str(
            'Extracted from "' + new_record.id + '" by easy_dna.extract_features()'
        )

        new_record.id = compute_id(seq_feature)
        new_record.name = new_record.id
        records = records + [new_record]

    return records


def write_records(key, records_dict, add_prefix=True):
    """Write a list of SeqRecords into Genbank files.
    """
    records = records_dict[key]
    root = flametree.file_tree(key)
    for j, record in enumerate(records):

        if add_prefix:
            filename_prefix = "feature_%s_" % j
        else:
            filename_prefix = ""
        record_name_alnum = "".join(x if x.isalnum() else "_" for x in record.name)
        record_filename = filename_prefix + record_name_alnum + ".gb"

        try:
            write_record(record, root._file(record_filename).open("w"), fmt="genbank")

        except Exception as err:
            print("Error writing", record_filename, str(err))


def make_part_dict(records_dict, min_sequence_length=20):
    """Make a full part list and a report.

    Uses records_dict by extractor_from_file() or extractor_from_batch().

    Parameters
    ==========
    records_dict
      Dictionary of sequence name: list of features as SeqRecords. 

    min_sequence_length
      Discard sequences with length less than this integer.
    """
    report_index = [
        "input_construct",
        "input_sequence",
        "sequence_string",
        "shared_with",
        "equal_to",
        "has_copy",
        "too_short",
        "notes",
    ]

    report = pd.DataFrame(
        data=None, index=report_index, columns=None, dtype=str, copy=False
    )
    common_parts_dict = dict()

    # This part is complex because it does two things, and will be simplified.
    # It makes a dictionary of all parts and a dataframe of part properties.
    # It also checks for a number of constraints and cases: sequence length,
    # shared sequences, shared names.
    for i, (key, record_list) in enumerate(records_dict.items()):
        for j, record in enumerate(record_list):
            s = pd.Series(None, index=report_index, name=i)
            s["input_construct"] = key
            s["input_sequence"] = record.name
            s["has_copy"] = False  # default

            if len(record) < min_sequence_length:
                s["too_short"] = True
                series_name = str(i) + "_" + str(j)
                report[series_name] = s
                continue
            else:
                s["too_short"] = False

            sequence_as_key = str(record.seq.lower())
            s["sequence_string"] = sequence_as_key

            if sequence_as_key in common_parts_dict.keys():
                s["has_copy"] = True
            else:
                remove = 1
                while record.name in report.loc["input_sequence"].array:
                    index = len(record.name) - remove
                    record.name = record.name[:index] + "X" + record.name[index + 1 :]
                    remove += 1

                    s["note"] = "renamed"

                s["input_sequence"] = record.name  # update name in record
                common_parts_dict[sequence_as_key] = record

            series_name = str(i) + "_" + str(j)
            report[series_name] = s

    report = report.T

    return (common_parts_dict, report)


def process_report(report):
    """Format the report prepared by make_part_dict().

    Finds common parts within constructs and identical sequences between
    constructs.
    """
    all_shared_with = pd.Series()
    all_equal_to = pd.Series()
    for sequence in report.loc[report["has_copy"]]["sequence_string"]:
        constructs = report.loc[report["sequence_string"] == sequence][
            "input_construct"
        ]
        shared_with_text = "_|_".join(constructs.unique())
        shared_with = pd.Series(
            data=shared_with_text, index=constructs.index, dtype=str
        )

        all_shared_with = pd.concat([all_shared_with, shared_with])

        parts = report.loc[report["sequence_string"] == sequence]["input_sequence"]
        equal_to_text = "_|_".join(parts.unique())
        equal_to = pd.Series(data=equal_to_text, index=parts.index, dtype=str)

        all_equal_to = pd.concat([all_equal_to, equal_to])

    report.loc[all_equal_to.index, "equal_to"] = all_equal_to.values
    report.loc[all_shared_with.index, "shared_with"] = all_shared_with.values

    return report
