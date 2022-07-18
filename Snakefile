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
configpath = Path(configfile)
cfg = Configuration(filepath=configpath)
global_params = cfg.global_params






rp_bwa_index = cfg.rule_params.BWA_Index_Reference
rule BWA_Index_Reference:
    input: global_params.reference_fasta
    output: 
        index_dir = global_params.files.reference_fasta_index_dir,
        valid_dir = "rc.BWA_Index_Reference.out"
    resources:
        **rp_bwa_index.resources,
        job_id=""
    envmodules:
        "gcc/6.2.0",
        "bwa/0.7.17"
    shell:
        "bwa index -p {output.index} {input.fasta_path}"
    

rule Align_FASTQ:
    input:
        fq1= resolve_fastq(config, keys="run={run}",fq1)
        fq2= resolve_fastq(config, keys="run={run}",fq2)
    shell:
        
