FROM python:3.7

RUN pip3 install synapseclient pandas numpy scipy sklearn

COPY validate_and_score.py /usr/local/bin/validate_and_score.py
COPY validate_sc1.py /usr/local/bin/validate_sc1.py
COPY score_sc1.py /usr/local/bin/score_sc1.py
COPY validate_sc2.py /usr/local/bin/validate_sc2.py
COPY score_sc2.py /usr/local/bin/score_sc2.py
