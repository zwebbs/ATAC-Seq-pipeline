# File Name: Snakefile
# Created By: ZW
# Created On: 2022-07-12
# Purpose: Runs the snakemake controller for the ATAC-Seq pipeline

# Module Imports
# ----------------------------------------------------------------------------

from gardnersnake import Configuration
from pathlib import Path


# Global Configuration and global outputs
# ----------------------------------------------------------------------------
configpath = Path(configfile)
cfg = Configuration(filepath=configpath)
global_pars = cfg.global_params


rule all:
    input:
        "rc.bwa_index_dir.out"

# Rule 0. Create BWA Index from FASTA Reference if not already present
# ----------------------------------------------------------------------------
rulepars_verify_index = cfg.rule_params.Verify_Index_Contents
rulepars_bwa_index = cfg.rule_params.BWA_Index_Reference
index_prefix = Path(global_params.files.reference_fasta_index_prefix).resolve()
index_dir = index_prefix.parent
index_files = [str(index_prefix / ext) for ext in (".amb",".ann",".bwt",".pac",".sa")]

# if the index directory already exists, check that the files are all there
# throw an error if not. I.e. if the index directory exists, it should be well
# formed and contain the necessary files. (the easiest way to achieve this is)
# to have a new subdirectory explicitly for the index
if index_dir.exists():
    rule Verify_Index_Contents:
    	input: manifest=index_files, checkdir=str(index_dir)
	output: rc="rc.bwa_index_dir.out"
    	resources: **rulepars_verify_index.resources, job_id=""
	shell:
	    "check_directory -o {output.rc} {input.manifest} {input.checkdir}"
else:
    rule BWA_Index_Reference:
	input: fasta_path = global_params.reference_fasta
        output: rc="rc.bwa_index_dir.out"
	params: index_dir=index_dir
    	resources: **rulepars_bwa_index.resources, job_id=""
    	envmodules: "gcc/6.2.0", "bwa/0.7.17"
        shell:
            "bwa index -p {params.index_dir} {input.fasta_path}"
    

