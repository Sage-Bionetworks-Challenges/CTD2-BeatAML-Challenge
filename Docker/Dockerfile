FROM python:3.7

RUN pip3 install synapseclient pandas numpy scipy sklearn lifelines

COPY sc*_utils.py /usr/local/bin/
COPY validate_sc*.py /usr/local/bin/
COPY score_sc*.py /usr/local/bin/
