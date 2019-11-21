#!/usr/bin/env python3

"""Validate SC2."""

import argparse
import json

import pandas as pd

SC2_DTYPE = {
    'lab_id': str,
    'responder': int
}


def main(submission, goldstandard, results):
    """Validate submission and write results to JSON.

    Args:
        submission: input file
        entity: Synapse entity type
        results: output file
    """

    invalid_reasons = []

    truth = pd.read_csv(goldstandard, dtype=SC2_DTYPE)
    try:
        predictions = pd.read_csv(submission, dtype=SC2_DTYPE)
    except Exception:
        invalid_reasons.append(
            "Prediction file must be a two-column CSV file.")

    else:

        # Column names are [lab_id, responder]
        header = predictions.columns.tolist()
        expected_cols = ["lab_id", "responder"]
        if header != expected_cols:
            invalid_reasons.append(
                f"Invalid columns names: {header}. Expected: {expected_cols}")

        else:

            # All expected lab_ids are accounted for
            pred_lab_ids = predictions.lab_id.tolist()
            gs_lab_ids = truth.lab_id.tolist()
            all_ids_found = all(
                lab_id in pred_lab_ids for lab_id in gs_lab_ids)
            if not all_ids_found:
                invalid_reasons.append(
                    "One or more expected lab_id(s) missing")

            # Responder values are Boolean (0, 1)
            responders = set(predictions.responder)
            if responders != {0, 1}:
                invalid_reasons.append(
                    "'responder' values should be a Boolean (0 or 1)")

    prediction_file_status = "INVALID" if invalid_reasons else "VALIDATED"

    result_dict = {'prediction_file_errors': "\n".join(invalid_reasons)[:500],
                   'prediction_file_status': prediction_file_status,
                   'round': 1}

    with open(results, 'w') as out:
        out.write(json.dumps(result_dict))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--submission_file",
                        required=True, help="Submission file")
    parser.add_argument("-g", "--goldstandard",
                        required=True, help="Truth file")
    parser.add_argument("-r", "--results",
                        required=True, help="Results file")

    args = parser.parse_args()
    main(args.submission_file, args.goldstandard, args.results)
