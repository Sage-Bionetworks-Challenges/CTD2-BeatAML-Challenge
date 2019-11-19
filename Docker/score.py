"""Implements validation and scoring functions for subchallenges 1 and 2."""

import numpy
import pandas


def listsAreEqual(l1, l2):
  """Returns true if the two lists have the same contents."""
  if len(l1) is not len(l2):
    return False
  return all([e1 == e2 for e1, e2 in zip(l1, l2)])


def validateSC1(submission_path, goldstandard_path):
  """Validate the format of SC1 submission file given the golden file.

  Args:
    submission_path: path to subchallenge 1 submission file
    goldstandard_path: path to the golden / ground truth file for subchallenge
        1.
  
  Raises:
    ValueError if submission file is invalid.
  """
  with open(submission_path, 'r') as s_file:
    with open(goldstandard_path, 'r') as g_file:
      validateSC1(s_file, g_file)


def validateSC1WithFiles(submission_file, goldstandard_file):
  """Same as above, but takes file handles."""
  try:
    submission = pandas.read_csv(submission_file)
  except Exception as e:
    raise ValueError("Not a properly formatted CSV") from e

  expected_columns = ['lab_id', 'drug', 'auc']
  if not listsAreEqual(expected_columns, submission.columns.tolist()):
    raise ValueError("Invalid columns. Got\n\t%s\nExpected\n\t%s" %
        (str(submission.columns.tolist()), str(expected_columns)))

  golden = pandas.read_csv(goldstandard_file)
  assert listsAreEqual(expected_columns, golden.columns.tolist()), (
      "Goldenfile has wrong columns.")


def scoreSC1(submission_path, goldstandard_path):
  """Returns primary and secondary scores for a SC1 submission.

  Args:
    submission_path: path to SC1 submission
    goldstandard_path: path to SC1 ground truth file.

  Raises:
    ValueError if submission file is invalid.

  Returns:
    A tuple with (spearman, pearson) correlation coefficients, the primary
    and secondary scores.
  """
  with open(submission_path, 'r') as s_file:
    with open(goldstandard_path, 'r') as g_file:
      scoreSC1WithFiles(s_file, g_file)


def scoreSC1WithFiles(submission_file, goldstandard_file):
  """Same as above, but takes file objects. Assumes files are validated."""
  dtype = {'lab_id': str, 'drug': str}
  submission = pandas.read_csv(submission_file, dtype=dtype)
  goldstandard = pandas.read_csv(goldstandard_file, dtype=dtype)

  # Create a single DataFrame which is indexed by (lab_id, drug) and has two
  # columns: [auc_submission, auc_goldstandard]
  indices = ['drug', 'lab_id']
  submission.set_index(indices, inplace=True)
  goldstandard.set_index(indices, inplace=True)
  joined = submission.join(
      goldstandard,
      lsuffix='_submission',
      rsuffix='_goldstandard')[['auc_submission', 'auc_goldstandard']]

  # Compute correlation for each drug.
  correlations = []
  for drug in joined.index.get_level_values('drug').unique():
    subset = joined.loc[drug]
    corr = subset.corr(method='spearman').auc_submission.auc_goldstandard
    correlations.append(corr)

  return numpy.mean(correlations)
