"""Test the validation and scoring functions for subchallenge 1.

Can be run with, eg:
  python3 -m unittest discover -f test/
"""

import pandas
import tempfile
import unittest

from score import scoreSC1WithFiles
from score import validateSC1WithFiles

class Subchallenge1Test(unittest.TestCase):
  def setUp(self):
    self.submission_df = pandas.DataFrame({
      'lab_id': ['1', '1', '2', '2', '3'],
      'drug': ['d1', 'd2', 'd1', 'd2', 'd1'],
      'auc': [0.1, 0.1, 0.2, 0.2, 0.3]
    })
    self.golden_df = pandas.DataFrame({
      'lab_id': ['1', '1', '2', '2', '3'],
      'drug': ['d1', 'd2', 'd1', 'd2', 'd1'],
      'auc': [0.1, 0.1, 0.2, 0.2, 0.3]
    })


  def runValidation(self):
    """Writes both dataframes to file and calls validateSC1."""
    with tempfile.TemporaryFile('w+') as submission_file:
      with tempfile.TemporaryFile('w+') as golden_file:
        self.submission_df.to_csv(submission_file, index=False)
        self.golden_df.to_csv(golden_file, index=False)
        submission_file.seek(0)
        golden_file.seek(0)

        validateSC1WithFiles(submission_file, golden_file)


  def runScoring(self):
    """Writes both dataframes to file and calls scoreSC1."""
    with tempfile.TemporaryFile('w+') as submission_file:
      with tempfile.TemporaryFile('w+') as golden_file:
        self.submission_df.to_csv(submission_file, index=False)
        self.golden_df.to_csv(golden_file, index=False)
        submission_file.seek(0)
        golden_file.seek(0)

        return scoreSC1WithFiles(submission_file, golden_file)


  def testSuccessfulValidation(self):
    self.runValidation()


  def testSuccessfulScore(self):
    self.runValidation()
    self.assertEqual(self.runScoring(), 1.0)

    self.submission_df.auc *= -1
    self.assertEqual(self.runScoring(), -1.0)

