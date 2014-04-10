#!/bin/sh

export PYTHON_VERSION=python_version
export PYTHONPATH=/pkg/lib/python${PYTHON_VERSION}/site-packages
export PYTHONPATH=${PYTHONPATH}:/opt/local/Library/Frameworks/Python.framework/Versions/${PYTHON_VERSION}lib/python${PYTHON_VERSION}/site-packages

scons -Q --tree=derived $@

