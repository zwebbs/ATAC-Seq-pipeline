# ATAC-Seq Pipeline


## Table of Contents


[Setting up the compute environment](#setting_up_compute)

[Running the pipeline](#run_pipeline)


\
<a name=setting_up_compute />
## Setting up the compute environment

To build the python environment for running the snakemake pipeline, load `python 3.8.10` or compatible equivalent, then create the virtual environment and install pip dependencies.
```bash
# load base python3
module load gcc/9.4.0   # only relavent for use on Gardner
module load python/3.8.10   # only relavent for use on Gardner

# create the virtual environment
python3 -m venv env/snakemake/

# activate the virtual envirionment and install package dependencies
source env/snakemake/bin/activate
python3 -m pip install -r env/snakemake/requirements.txt

```

\
<a name=run_pipeline />
## Running the pipeline

The workflow has a controller script named `run_pipeline.sh` that can be used to run the pipeline on Gardner or a similar HPC infastruture. There are three argument flags for the controller script: `-j` which specifies the number of concurrent jobs allowed to be submitted to the cluster, `-c` which specifies the configuration file for the workflow, and (optionally) `-d` which indicates that the snakemake workflow should be run in dry run mode (equivalent to using the snakemake `--dry-run` flag).

```bash
# run the snakemake pipeline, optionally as a dry run --using the -d flag
./run_pipeline.sh -j 24 -c config/template-config.yaml # -d
```



