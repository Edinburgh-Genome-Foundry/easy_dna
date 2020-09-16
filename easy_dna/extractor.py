import pandas as pd
import flametree
from .io import load_record, records_from_data_files, write_record


def extract_from_input(
    filename=None,
    directory=None,
    construct_list=None,
    direct_sense=True,
    output_path=None,
    min_sequence_length=20,
):
    """Extract features from input and return in a dictionary.

    Optionally save the features in separate files.

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

    output_path
      Path for the exported feature and report files.

    min_sequence_length
      Discard sequences with length less than this integer.
    """
    genbank_id_limit = 20  # GenBank format hard limit for name
    if construct_list:
        pass
    elif filename:
        input_record = load_record(
            filename, record_id="auto", upperize=False, id_cutoff=genbank_id_limit
        )
        construct_list = [input_record]
    elif directory:
        construct_list = records_from_data_files(filepaths=None, folder=directory)
    else:
        raise TypeError("Specify one of 'construct_list', 'filename' or 'directory'.")

    records_dict = dict()
    recordname_list = []
    for input_record in construct_list:
        records = extract_features(input_record, direct_sense)
        record_name = input_record.name[0:genbank_id_limit]
        # This part makes the key (used as dir name) unique by appending a copynumber:
        number_of_name_occurrences = recordname_list.count(record_name)
        if number_of_name_occurrences:
            key = "%s_%s" % (record_name, number_of_name_occurrences + 1)
        else:
            key = record_name

        recordname_list.append(record_name)

        records_dict[key] = records

    parts_report = make_part_dict(records_dict, min_sequence_length=min_sequence_length)
    processed_report = process_report(parts_report[1])

    all_parts_dict = parts_report[0]
    records_dict["all_parts"] = list(all_parts_dict.values())

    if output_path is not None:
        root = flametree.file_tree(output_path)

        for key, records in records_dict.items():

            record_dir = root._dir(key)

            record_name_alnum_list = []
            for record in records:

                record_name_alnum = "".join(
                    x if x.isalnum() else "_" for x in record.name
                )
                # This part makes the filename unique by appending a copynumber:
                number_of_occurrences = record_name_alnum_list.count(record_name_alnum)
                if number_of_occurrences:
                    record_filename = "%s_%s.gb" % (
                        record_name_alnum,
                        number_of_occurrences + 1,
                    )
                else:
                    record_filename = record_name_alnum + ".gb"

                record_name_alnum_list.append(record_name_alnum)

                record_file_path = record_dir._file(record_filename)

                try:
                    write_record(record, record_file_path, fmt="genbank")

                except Exception as err:
                    print("Error writing", record_filename, str(err))

        processed_report.to_csv(root._file("report.csv").open("w"))

    records_dict["processed_report"] = processed_report

    return records_dict


def compute_id(seq_feature):
    """Computes an ID for extract_features().

    The ID can be used as a SeqRecord/Genbank name or id attribute.
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
        for other_seq in (
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


def make_part_dict(records_dict, min_sequence_length=20):
    """Make a full part list and a report.

    Uses records_dict made by run_extraction().

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
        "shared_with",
        "equal_to",
        "has_copy",
        "too_short",
        "notes",
        "sequence_string",
    ]

    report = pd.DataFrame(
        data=None, index=report_index, columns=None, dtype=str, copy=False
    )
    all_parts_dict = dict()

    # This part is complex because it does two things, and could be simplified.
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

            if sequence_as_key in all_parts_dict.keys():
                s["has_copy"] = True
            else:
                remove = 1
                while record.name in report.loc["input_sequence"].array:
                    index = len(record.name) - remove
                    record.name = record.name[:index] + "X" + record.name[index + 1 :]
                    remove += 1

                    s["note"] = "renamed"
                record.id = record.name

                s["input_sequence"] = record.name  # update name in record
                all_parts_dict[sequence_as_key] = record

            series_name = str(i) + "_" + str(j)
            report[series_name] = s

    report = report.T

    return (all_parts_dict, report)


def process_report(report):
    """Format the report prepared by make_part_dict().

    The function finds common parts within constructs and identical sequences
    between constructs.
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

    report.drop("has_copy", axis=1, inplace=True)

    return report
