#!/bin/sh
cd `dirname $0`
export PATH=$PATH:/home/jab69/virtualenvs/pysb/bin
git pull --ff-only && python submit_build.py
