#!/bin/bash

CONDA_ENV_PY3=chip-nexus-pipeline
CONDA_ENV_PY2=chip-nexus-pipeline-python2
CONDA_ENV_OLD_PY3=chip-nexus-pipeline-python3

conda env remove -n ${CONDA_ENV_PY3} -y
conda env remove -n ${CONDA_ENV_PY2} -y
conda env remove -n ${CONDA_ENV_OLD_PY3} -y

