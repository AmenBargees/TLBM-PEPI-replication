#!/bin/csh -f

set file = "/home/peter/src/tlbm/projects/`cat /home/peter/src/tlbm/PROJECT`/`cat /home/peter/src/tlbm/RUN`.in"
grep $1 $file
