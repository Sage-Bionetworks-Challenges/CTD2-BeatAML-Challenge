#!/usr/bin/env python3

"""Validate SC2."""

import argparse
import json

import pandas as pd


def main(submission, entity_type, results="results.json"):
    """Validate submission and write results to JSON.

    Args:
        submission: input file
        entity: Synapse entity type
        results: output file
    """

    invalid_reasons = []

    if not submission.endswith(".csv"):
        invalid_reasons.append(
            "Prediction file must be a two-column CSV file.")
    else:
        predictions = pd.read_csv(submission, sep=",")

        header = list(predictions.columns)
        try:
            assert header == ["lab_id", "responder"]
        except AssertionError:
            invalid_reasons.append("Column names must be: lab_id,responder.")

        responses = set(predictions.responder)
        try:
            assert responses == {0, 1}
        except AssertionError:
            invalid_reasons.append(
                " Responder values must be a Boolean (0 or 1).")

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
    parser.add_argument("-e", "--entity_type",
                        required=True, help="Synapse entity type")
    parser.add_argument("-r", "--results",
                        required=True, help="Results file")

    args = parser.parse_args()
    main(args.submission_file)
