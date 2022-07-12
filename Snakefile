# File Name: Snakefile
# Created By: ZW
# Created On: 2022-07-12
# Purpose: Runs the snakemake controller for the ATAC-Seq pipeline
# 


# rule Global:

rpars_bwa_index
rule BWA_Index_Reference:
    input: fasta_path = resolve_file(config, keys="", "reference_fasta",exists=True)
    output: 
        index = resolve_file(config, keys="", referenceexists=False)    
        valid_dir = "rc.BWA_Index_Reference.out"
    
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
        
