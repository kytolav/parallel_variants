#!/bin/bash
export FLASK_APP=app/run_tracker/run_tracker/__init__.py
export FLASK_DEBUG=1
python3 -m flask run --port $1