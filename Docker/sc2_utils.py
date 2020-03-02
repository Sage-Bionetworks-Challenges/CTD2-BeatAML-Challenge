"""Utilities for validating and scoring SC2."""

import lifelines
import pandas
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv


SC2_SUBMISSION_DTYPE = {
    'lab_id': str,
    'survival': 'float64'
}


SC2_GOLDEN_DTYPE = {
    'lab_id': str,
    'vitalStatus': str,
    'overallSurvival': 'float64'
}


def validateSC2(submission, goldstandard):
    """Returns a results dict."""
    invalid_reasons = []

    truth = pandas.read_csv(goldstandard, dtype=SC2_GOLDEN_DTYPE)
    try:
        predictions = pandas.read_csv(submission, dtype=SC2_SUBMISSION_DTYPE)
    except Exception as e:
        invalid_reasons.append(
            "Prediction file is improperly formatted. Must have two columns, %s. Error: %s"
            % (str(SC2_SUBMISSION_DTYPE), str(e)))

    else:
        header = predictions.columns.tolist()
        expected_cols = ["lab_id", "survival"]
        if header != expected_cols:
            invalid_reasons.append(
                f"Invalid column names: {header}. Expected: {expected_cols}")

        else:
            # All ground truth lab_ids have predictions.
            pred_lab_ids = set(predictions.lab_id)
            gs_lab_ids = set(truth.lab_id)
            missing_ids = gs_lab_ids - pred_lab_ids

            if missing_ids:
                invalid_reasons.append(
                    f"Missing {len(missing_ids)} predictions.")

            if predictions.lab_id.duplicated().any():
                invalid_reasons.append(
                    "There are duplicate predictions.")

            # Survival should be... any non null number.
            num_nan_survivals = predictions.survival.isna().sum()
            if num_nan_survivals:
                invalid_reasons.append(
                    f"{num_nan_survivals} predictions were NAN.")

    prediction_file_status = "INVALID" if invalid_reasons else "VALIDATED"

    result_dict = {'prediction_file_errors': "\n".join(invalid_reasons)[:500],
                   'prediction_file_status': prediction_file_status,
                   'round': 1}
    return result_dict


def responseToSurvivalMatrix(response):
    """Converts a response.csv to a survival matrix expected by scikitsurv."""
    return Surv.from_dataframe('vitalStatus', 'overallSurvival',
        pandas.concat([
            (response.vitalStatus == 'Dead'),
            response.overallSurvival
        ], axis=1))


def scoreSC2(submission, goldstandard, trainingdata):
    """Score a subchallenge 2 submission, returns concordance index."""
    predictions = pandas.read_csv(submission, dtype=SC2_SUBMISSION_DTYPE)

    truth = pandas.read_csv(goldstandard, dtype=SC2_GOLDEN_DTYPE)

    joined = (truth
        .set_index('lab_id')
        .join(
          predictions
          .set_index('lab_id')
          .rename(columns={'survival': 'prediction'})
        ))

    ci = lifelines.utils.concordance_index(
        joined.overallSurvival,
        joined.prediction,
        (joined.vitalStatus == 'Dead'))

    DAYS = 365
    trainingdata = pandas.read_csv(trainingdata)

    auc = cumulative_dynamic_auc(
        responseToSurvivalMatrix(trainingdata),
        responseToSurvivalMatrix(joined[['vitalStatus', 'overallSurvival']]),
        -joined.prediction.to_numpy(),
        [DAYS]
    )[0][0]

    return (ci, auc)
