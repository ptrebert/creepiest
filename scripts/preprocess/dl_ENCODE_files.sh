#!/bin/bash

TARGET_DIR=/TL/deep/fhgfs/projects/pebert/thesis/biodata/dlfolder/encode

LISTING=/TL/deep/fhgfs/projects/pebert/thesis/biodata/dlfolder/listing_encode.txt

cd ${TARGET_DIR}

CMD="xargs -n 1 -a ${LISTING} curl -L -O -s -S"

${CMD}

exit $?

