#!/bin/bash

# File Name: run_pipeline.sh
# Created On: 2022-07-07
# Created By: ZW
# Purpose: runs ATAC-Seq alignment, peak calling, and analysis
# pipeline for the given sets of FASTQ data.

# check passed commandline arguments
## of which there are three.
###   -j (Number of concurrent jobs) [required]
###   -c <analysis configuration file .yaml> [required]
###   -d (BOOLEAN flag to complete a snakemake dry run) [optional]

PIPELINE_NAME="EHT-Test-ATAC-Seq"

while getopts ":j:c:d" 'opt';
do
    case ${opt} in
        j) 
            parallel_jobs=${OPTARG}
            ;;  # capture command line argument for number of possible concurrent jobs
        c) 
            config_file=${OPTARG}
            ;;  # capture command line argument for YAML config file
        d) 
            dry_run_flag="--dry-run"
            ;;  # capture boolean for whether or not its a 'dry run'
        *)
            echo "INVALID_OPTION -- ${OPTARG}"
            exit 1
            ;;  # capture all other (invalid) inputs
    esac
done


# print the passed commandline options
echo "Number of concurrent jobs: ${parallel_jobs}"
echo "Selected config file: ${config_file}"
echo "dry run ?: ${dry_run_flag}"


# load modules
module load python/cpython-3.8.5
module load java/1.8
module load samtools/1.9
module load R/4.1.0


# set script globals
export PICARD="/home/zweber/tools/gatk-4.2.6.1/picard.jar"

# run the snakemake workflow
snakemake --snakefile Snakefile --use-envmodules \
    -j ${parallel_jobs} -kp --rerun-incomplete \
    --default-resources mem_mb=64000 disk_mb=100000 \
    --config yaml_config=${config_file} \
    --cluster "sbatch --job-name={rulename}_{resources.job_id} \
     --partition=broadwl \
     --error={resources.logs}/{rulename}_{resources.job_id}.err \
     --output={resources.logs}/{rulename}_{resources.job_id}.out \
     --nodes={resources.nodes} \
     --ntasks-per-node={resources.cpus_per_node} \
     --mem={resources.total_memory_mb} \
     --time={resources.walltime}" \
     ${dry_run_flag} 1>  "${PIPELINE_NAME}.out" 2> "${PIPELINE_NAME}.err"

# write that the pipeline is complete
echo "-----------------------------"
echo "      ${PIPELINE_NAME}       "
echo "      Pipeline Done!         "

