#!/usr/bin/env python3

"""Validate SC1."""

import argparse
import json

import validate_and_score


def main(submission, goldstandard, results):
    """Validate submission and write results to JSON.

    Args:
        submission: input file
        entity: Synapse entity type
        results: output file
    """

    invalid_reasons = ""

    try:
        validate_and_score.validateSC1(submission, goldstandard)
    except (AssertionError, ValueError) as e:
        invalid_reasons = str(e)

    prediction_file_status = "INVALID" if invalid_reasons else "VALIDATED"

    result_dict = {'prediction_file_errors': invalid_reasons[:500],
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
