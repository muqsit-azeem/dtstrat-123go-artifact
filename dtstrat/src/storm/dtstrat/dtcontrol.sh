#!/bin/bash
#
# Run the following steps in storm's working directory:
#   python -m venv venv
#   source venv/bin/activate
#   pip install dtcontrol
#

. venv/bin/activate
dtcontrol "$@"
