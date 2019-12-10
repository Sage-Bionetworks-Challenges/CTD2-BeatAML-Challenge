"""Implements validation and scoring functions for subchallenge 1."""

import numpy
import pandas


SC1_DTYPE = {
  'lab_id': str,
  'inhibitor': str,
  'auc': numpy.float64
}


def setsAreEqual(l1, l2):
  """Returns true if the two lists have the same contents."""
  if len(l1) is not len(l2):
    return False
  return all([e1 == e2 for e1, e2 in zip(sorted(l1), sorted(l2))])


def validateSC1(submission_path, goldstandard_path):
  """Validate the format of SC1 submission file given the golden file.

  Args:
    submission_path: path to subchallenge 1 submission file
    goldstandard_path: path to the golden / ground truth file for subchallenge
        1.
  
  Raises:
    ValueError if submission file is invalid.
  """
  try:
    submission = pandas.read_csv(submission_path, dtype=SC1_DTYPE)
  except Exception as e:
    raise ValueError(
      "Prediction file is improperly formatted. Must have two columns, %s. Error: %s"
      % (str(SC1_DTYPE), str(e)))

  expected_columns = ['lab_id', 'inhibitor', 'auc']
  if not setsAreEqual(expected_columns, submission.columns.tolist()):
    raise ValueError("Invalid columns. Got\n\t%s\nExpected\n\t%s" %
        (str(submission.columns.tolist()), str(expected_columns)))

  golden = pandas.read_csv(goldstandard_path, dtype=SC1_DTYPE)
  assert setsAreEqual(expected_columns, golden.columns.tolist()), (
      "Goldenfile has wrong columns.")

  # Check for duplicates.
  duplicates = submission.duplicated(subset=['lab_id', 'inhibitor'])
  if duplicates.any():
    raise ValueError("Found %d duplicate(s) including for (%s,%s)." % (
      duplicates.sum(),
      submission[duplicates].lab_id.iloc[0],
      submission[duplicates].inhibitor.iloc[0]))

  # Check for 1:1 correspondence between entries, by computing the set
  # difference between two multiindices.
  indices = ['inhibitor', 'lab_id']
  submission.set_index(indices, inplace=True)
  golden.set_index(indices, inplace=True)
  extra_rows = submission.index.difference(golden.index)
  if not extra_rows.empty:
    raise ValueError("Found %d unexpected row(s) in submission." %
      extra_rows.shape[0])

  missing_rows = golden.index.difference(submission.index)
  if not missing_rows.empty:
    raise ValueError("Missing %d row(s) in submission." %
      missing_rows.shape[0])


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
  submission = pandas.read_csv(submission_path, dtype=SC1_DTYPE)
  goldstandard = pandas.read_csv(goldstandard_path, dtype=SC1_DTYPE)

  # Create a single DataFrame which is indexed by (lab_id, inhibitor) and has two
  # columns: [auc_submission, auc_goldstandard]. For both submission and 
  # goldstandard, we create a multindex with two dimensions: inhibitor and lab_id,
  # then do a join-by-index.
  indices = ['inhibitor', 'lab_id']
  submission.set_index(indices, inplace=True)
  goldstandard.set_index(indices, inplace=True)
  joined = submission.join(
      goldstandard,
      lsuffix='_submission',
      rsuffix='_goldstandard')[['auc_submission', 'auc_goldstandard']]

  # Compute correlation for each inhibitor.
  spearmans = []
  pearsons = []
  for inhibitor in joined.index.get_level_values('inhibitor').unique():
    subset = joined.loc[inhibitor]
    spearmans.append(
        subset.corr(method='spearman').auc_submission.auc_goldstandard)
    pearsons.append(
        subset.corr(method='pearson').auc_submission.auc_goldstandard)

  return (numpy.mean(spearmans), numpy.mean(pearsons))
