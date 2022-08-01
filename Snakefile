# File Name: Snakefile
# Created By: ZW
# Created On: 2022-07-12
# Purpose: Runs the snakemake controller for the ATAC-Seq pipeline

# Module Imports
# ----------------------------------------------------------------------------
from gardnersnake import Configuration
from pathlib import Path


# Global Configuration
# ----------------------------------------------------------------------------
yaml_config_filepath = Path(config["yaml_config"])
cfg = Configuration(filepath=yaml_config_filepath)
cfg.load()

GLOBALS = cfg.global_params
WORKDIR = GLOBALS.working_directory
REF = Path(GLOBALS.files.reference_fasta)

# set analysis workdir
workdir: WORKDIR

# build a map of the sequencing data
# TODO: implement map that suits ATAC data
#


# Rule 0. Pipeline Global Returns
# ----------------------------------------------------------------------------
rule All:
    input: (REF.parent / "BWA_Index_Reference.rc.out")


# Rule 1. Create BWA Index from FASTA Reference
# -----------------------------------------------------------------------------
bwa_index_rp = cfg.get_rule_params(rulename="BWA_Index_Reference")
rule BWA_Index_Reference:
    input: fasta_path = REF
    params:  **(bwa_index_rp.parameters), index_dir = REF.parent
    resources:  **(bwa_index_rp.resources), job_id = "glob"
    envmodules: *(bwa_index_rp.parameters.envmodules)
    output: rc = (REF.parent / "BWA_Index_Reference.rc.out")
    shell:
        "bwa index {input.fasta_path}"
        " && check_directory -o {output.rc}"
        " {params.index_files} {params.index_dir}"

