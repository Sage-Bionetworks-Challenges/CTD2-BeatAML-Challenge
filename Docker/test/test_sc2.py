"""Test the validation and scoring functions for subchallenge 2.

Can be run with, eg:
  python3 -m unittest discover -f test/
"""

import json
import tempfile
import unittest

import numpy
import pandas

from sc2_utils import validateSC2
from sc2_utils import scoreSC2


class Subchallenge2Test(unittest.TestCase):
  def setUp(self):
    self.submission_df = pandas.DataFrame({
      'lab_id': ['1', '2', '3'],
      'survival': [0, 2.2, 3.3]
    })
    self.golden_df = pandas.DataFrame({
      'lab_id': ['2', '3'],
      'vitalStatus': ['Dead', 'Dead'],
      'overallSurvival': [300, 400],
    })
    self.trainingdata_df = pandas.DataFrame({
      'lab_id': ['t1', 't2', 't3'],
      'vitalStatus': ['Dead', 'Dead', 'Dead'],
      'overallSurvival': [200, 300, 400],
    })
    self.tmp_dir = tempfile.TemporaryDirectory()


  def tearDown(self):
    self.tmp_dir.cleanup()


  def runValidation(self, expected_error_str=''):
    """Writes both dataframes to file and calls validateSC1."""
    submission_fname = str(self.tmp_dir.name) + '/submission.csv'
    golden_fname = str(self.tmp_dir.name) + '/golden.csv'
    self.submission_df.to_csv(submission_fname, index=False)
    self.golden_df.to_csv(golden_fname, index=False)

    results = validateSC2(submission_fname, golden_fname)
    if not expected_error_str:
      self.assertEqual('', results['prediction_file_errors'])
    else:
      self.assertEqual('INVALID', results['prediction_file_status'])
      self.assertRegex(results['prediction_file_errors'], expected_error_str)


  def runScoring(self):
    """Writes both dataframes to file and calls scoreSC1."""
    submission_fname = str(self.tmp_dir.name) + '/submission.csv'
    golden_fname = str(self.tmp_dir.name) + '/golden.csv'
    trainingdata_fname = str(self.tmp_dir.name) + '/trainingdata.csv'
    self.submission_df.to_csv(submission_fname, index=False)
    self.golden_df.to_csv(golden_fname, index=False)
    self.trainingdata_df.to_csv(trainingdata_fname, index=False)

    return scoreSC2(submission_fname, golden_fname, trainingdata_fname)


  def testSuccessfulValidatation(self):
    self.runValidation()


  def testValidateMissingColumns(self):
    self.submission_df.drop(columns=['survival'], inplace=True)
    self.runValidation(expected_error_str='Invalid column names')


  def testValidateMissingSpecimen(self):
    self.submission_df.drop(1, inplace=True)
    self.runValidation(expected_error_str='Missing 1 prediction')


  def testValidateNanValue(self):
    self.submission_df = self.submission_df.append(
        {'lab_id': '4', 'survival': numpy.NAN},
        ignore_index=True)
    self.runValidation(expected_error_str='predictions were NAN')


  def testValidateDuplicatePrediction(self):
    self.submission_df = self.submission_df.append(
        {'lab_id': '3', 'survival': numpy.NAN},
        ignore_index=True)
    self.runValidation(expected_error_str='There are duplicate predictions.')


  def testValidateMisplacedString(self):
    self.submission_df = self.submission_df.append(
        {'lab_id': '3', 'survival': 'I AM A BAD STRING'},
        ignore_index=True)
    self.runValidation(expected_error_str='could not convert string to float')


  def testSuccessfulScoring(self):
    (ci, auc) = self.runScoring()
    self.assertEqual(1.0, ci)
    self.assertEqual(1.0, auc)


  def testAnAlmostPerfectPrediction(self):
    """All pairs are correct except for the last one."""
    self.submission_df = pandas.DataFrame({
      'lab_id': [str(v) for v in range(5)],
      'survival': [1, 2, 3, 5, 4]  # <--- Note the one pair swap.
    })
    self.golden_df = pandas.DataFrame({
      'lab_id': [str(v) for v in range(5)],
      'vitalStatus': ['Dead'] * 5,
      'overallSurvival': [ii * 100 for ii in range(5)]
    })
    num_pairs = 5 * 4 / 2
    (ci, auc) = self.runScoring()
    self.assertEqual((num_pairs - 1.0) / num_pairs, ci)
    self.assertEqual(auc, 0.75)
