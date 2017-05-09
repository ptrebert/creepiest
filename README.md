# The CREEPIEST tool

The CRoss spEcies EPIgenome ESTimation tool

## Code maturity

Stable beta :neutral_face:

## Installation

At the moment, no install packages/routines are provided.

Create an environment as specified in `config/environment/conda_env_creepiest.yml` and
get current source from this repository. Set `PATH` and `PYTHONPATH` as appropriate
for your system and check if the CREEPIEST tool is executable via `creepiest.py --help`.

## Functional modules/commands

- apply: apply a trained model to test data
- compfeat: compute features for data
- convert-bedgraph: convert bedGraph file to HDF5
- convert-map: convert pairwise map (alignment) to index file in HDF5 format
- convert-region: convert BED-like file to HDF5
- dump: dump HDF5 to text-based file (works for index and region files)
- info: print basic HDF5 structure to screen
- map_signal: map signal tracks in HDF5 format from target to query
- match: for a given set of genomic regions, find similar regions in the genomic complement
- merge: fuse together datasets based on region names
- norm: perform quantile normalization
- train: train a machine learning model supporting a scikit-learn interface

## Dysfunctional modules/commands

- convert-chain: deprecated; will be removed at some point
- convert-motifdb: convert FIMO output to fixed-sized blocks for faster access; deprecated
- correlation: compute genome-wide correlation between signal tracks; does not yet use current index files, will
 be updated at some point
- tests: unit tests still use old index format, will be updated at some point