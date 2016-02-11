#!/bin/bash

TARGET_DIR=/TL/deep/fhgfs/projects/pebert/crppipe/encode/rawDL

LISTINGHG19=/TL/deep/fhgfs/projects/pebert/crppipe/encode/hg19_files.txt

LISTINGMM9=/TL/deep/fhgfs/projects/pebert/crppipe/encode/mm9_files.txt

cd ${TARGET_DIR}

CMD="xargs -n 1 -a ${LISTINGHG19} curl -L -O -s -S"

${CMD}

CMD="xargs -n 1 -a ${LISTINGMM9} curl -L -O -s -S"

${CMD}

exit $?

