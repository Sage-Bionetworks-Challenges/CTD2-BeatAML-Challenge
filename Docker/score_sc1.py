#!/usr/bin/env python3

"""Score SC1."""

import argparse
import json

import validate_and_score


def main(submission, goldstandard, results):
    """Get scores and write results to json

    Args:
        submissionfile: Participant submission file path
        goldstandard: Goldstandard file path
        results: File to write results to
        path_to_treecmp: Path to TreeCmp
    """
    score_dict = {
        'prediction_file_status': "SCORED"
    }

    spearman, pearson = validate_and_score.scoreSC1(submission, goldstandard)
    score_dict['spearman'] = spearman
    score_dict['pearson'] = pearson

    with open(results, 'w') as output:
        output.write(json.dumps(score_dict))


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
