#!/bin/sh
cd `dirname $0`
git pull --ff-only && bsub -q short -W 1:00 make clean html
