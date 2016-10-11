#!/bin/bash

INSTALLPATH="/TL/epigenetics2/work/pebert/conda/envs/creepiest"

cd ${INSTALLPATH} &&
mkdir -p ./etc/conda/activate.d &&
mkdir -p ./etc/conda/deactivate.d &&
touch ./etc/conda/activate.d/env_vars.sh &&
touch ./etc/conda/deactivate.d/env_vars.sh

CRPBINPATH="/home/pebert/work/code/mpggit/creepiest"
CRPPYTHONPATH="/home/pebert/work/code/mpggit/creepiest"

echo "#!/bin/sh" > ./etc/conda/activate.d/env_vars.sh

echo "" >> ./etc/conda/activate.d/env_vars.sh

echo "export PATH=${CRPBINPATH}:\$PATH" >> ./etc/conda/activate.d/env_vars.sh

echo "export PYTHONPATH=${CRPPYTHONPATH}:\$PYTHONPATH" >> ./etc/conda/activate.d/env_vars.sh



