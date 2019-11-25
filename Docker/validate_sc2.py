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

            # There is 1:1 mapping relationship between lab_ids in gs
            # and predictions (all expected ids found, no dups, etc.)
            pred_lab_ids = predictions.lab_id.tolist()
            uniq_pred_ids = set(pred_lab_ids)
            gs_lab_ids = set(truth.lab_id.tolist())

            missing_ids = gs_lab_ids - uniq_pred_ids
            unknown_ids = uniq_pred_ids - gs_lab_ids

            if missing_ids:
                invalid_reasons.append(
                    f"Missing lab_ids: {sorted(missing_ids)}")

            if unknown_ids:
                invalid_reasons.append(
                    f"Unknown lab_ids found: {sorted(unknown_ids)}")

            if len(pred_lab_ids) != len(uniq_pred_ids):
                invalid_reasons.append("Some lab_ids are duplicated")

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
