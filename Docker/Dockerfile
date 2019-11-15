FROM python:3

RUN pip install synapseclient pandas scipy sklearn

COPY validate.py /usr/local/bin/validate.py
COPY score.py /usr/local/bin/score.py
