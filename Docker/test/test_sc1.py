"""Test the validation and scoring functions for subchallenge 1.

Can be run with, eg:
  python3 -m unittest discover -f test/
"""

import pandas
import tempfile
import unittest

from validate_and_score import scoreSC1
from validate_and_score import validateSC1

class Subchallenge1Test(unittest.TestCase):
  def setUp(self):
    self.submission_df = pandas.DataFrame({
      'lab_id': ['1', '1', '2', '2', '3'],
      'inhibitor': ['d1', 'd2', 'd1', 'd2', 'd1'],
      'auc': [0.1, 0.1, 0.2, 0.2, 0.3]
    })
    self.golden_df = pandas.DataFrame({
      'lab_id': ['1', '1', '2', '2', '3'],
      'inhibitor': ['d1', 'd2', 'd1', 'd2', 'd1'],
      'auc': [0.1, 0.1, 0.2, 0.2, 0.3]
    })
    self.tmp_dir = tempfile.TemporaryDirectory()


  def tearDown(self):
    self.tmp_dir.cleanup()


  def runValidation(self):
    """Writes both dataframes to file and calls validateSC1."""
    submission_fname = str(self.tmp_dir.name) + '/submission.csv'
    golden_fname = str(self.tmp_dir.name) + '/golden.csv'
    self.submission_df.to_csv(submission_fname, index=False)
    self.golden_df.to_csv(golden_fname, index=False)

    validateSC1(submission_fname, golden_fname)


  def runScoring(self):
    """Writes both dataframes to file and calls scoreSC1."""
    submission_fname = str(self.tmp_dir.name) + '/submission.csv'
    golden_fname = str(self.tmp_dir.name) + '/golden.csv'
    self.submission_df.to_csv(submission_fname, index=False)
    self.golden_df.to_csv(golden_fname, index=False)

    return scoreSC1(submission_fname, golden_fname)


  def testSuccessfulValidation(self):
    self.runValidation()


  def testValidationCatchesMissingColumn(self):
    self.submission_df.drop(columns=['inhibitor'], inplace=True)
    with self.assertRaisesRegex(ValueError, 'Invalid columns') as cm:
      self.runValidation()


  def testValidationCatchesBadDtype(self):
    self.submission_df.auc = self.submission_df.astype('object')
    self.submission_df.auc.at[0] = 'bad_entry'
    with self.assertRaisesRegex(ValueError, 'Not a properly formatted CSV'):
      self.runValidation()


  def testValidationCatchesDuplicates(self):
    self.submission_df = self.submission_df.append(
        {'inhibitor': 'd1', 'lab_id': '1', 'auc': 0.12345},
        ignore_index=True)
    with self.assertRaisesRegex(ValueError, r'1 duplicate.*1.*d1'):
      self.runValidation()


  def testValidationCatchesExtraRow(self):
    # Use a new inhibitor, for no collisions.
    self.submission_df = self.submission_df.append(
        {'inhibitor': 'd4', 'lab_id': '1', 'auc': 0.12345},
        ignore_index=True)
    with self.assertRaisesRegex(ValueError, r'1 unexpected row'):
      self.runValidation()


  def testValidationCatchesMissingRow(self):
    self.submission_df = self.submission_df.drop([0])
    with self.assertRaisesRegex(ValueError, r'Missing 1 row'):
      self.runValidation()

  
  def testValidationAcceptsDifferentColumnOrder(self):
    self.submission_df = self.submission_df[['inhibitor', 'lab_id', 'auc']]
    self.runValidation()


  def testSuccessfulScore(self):
    self.runValidation()
    self.assertEqual(self.runScoring(), (1.0, 1.0))

    self.submission_df.auc *= -1
    self.assertEqual(self.runScoring(), (-1.0, -1.0))


  def testScoreSucceedsWithWeirdColumnOrdering(self):
    self.submission_df = self.submission_df[['inhibitor', 'lab_id', 'auc']]
    self.assertEqual(self.runScoring(), (1.0, 1.0))


