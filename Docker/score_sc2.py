#!/usr/bin/env python3
"""Score SC2."""

import argparse
import json

from sc2_utils import scoreSC2


def main(submission, goldstandard, trainingdata, results):
    """Get scores and write results to json

    Args:
        submissionfile: Participant submission file path
        goldstandard: Goldstandard file path
        trainingdata: Training data responses file path
        results: File to write results to
    """
    (ci, auc) = scoreSC2(submission, goldstandard, trainingdata)
    score_dict = {
        'prediction_file_status': "SCORED",
        'concordance_index': ci,
        'auc': auc
    }

    with open(results, 'w') as output:
        output.write(json.dumps(score_dict))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--submission_file",
                        required=True,
                        help="Submission file")
    parser.add_argument("-g",
                        "--goldstandard",
                        required=True,
                        help="Truth file")
    parser.add_argument("-t",
                        "--trainingdata",
                        required=True,
                        help="File with responses from the training data.")
    parser.add_argument("-r", "--results", required=True, help="Results file")

    args = parser.parse_args()
    main(args.submission_file,
        args.goldstandard,
        args.trainingdata,
        args.results)
