#!/bin/sh
cd `dirname $0`
git pull --ff-only && bsub -q short -W 6:00 make clean html
