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

# define global configuration objects
GLOBALS = cfg.global_params
WORKDIR = GLOBALS.working_directory
workdir: WORKDIR

# define directories and files objects
DIRECTORIES = GLOBALS.misc.directories
FILES = GLOBALS.misc.files

# define specific directories that are used often
SCRATCH_DIR = Path(DIRECTORIES.scratch_directory)
ALIGN_DIR = Path(DIRECTORIES.align_directory)
METRICS_DIR = Path(DIRECTORIES.metrics_directory)
SCRIPTS_DIR = Path(DIRECTORIES.scripts_directory)
PLOTS_DIR = Path(DIRECTORIES.plots_directory)
WIG_DIR = Path(DIRECTORIES.wig_directory)
PEAKS_DIR = Path(DIRECTORIES.peaks_directory)
LOGS_DIR = Path(DIRECTORIES.logs_directory)

# define specific files that are used often
REF = FILES.reference_fasta
SEQ = FILES.sequencing

# build a map of the sequencing data
RUN_IDS = [seq_run["run_id"] for seq_run in SEQ]
EXP_IDS = [seq_run["experiment_id"] for seq_run in SEQ]
EXPS = list(set(EXP_IDS))
SEQ_MAP = dict(zip(RUN_IDS, SEQ))

# set chromosomes for analysis
ALLOWED_CHRS = ["chr" + str(i) for i in [*list(range(1,23)),"X"]]

# inline utility function for getting fastq files
def get_fastq_file(run_id, pe_orient):
    return SEQ_MAP[run_id][f"reads_{pe_orient}"]


# Rule 0. Pipeline Global Returns
# ----------------------------------------------------------------------------
rule All:
    input: 
        (REF.parent / "BWA_Index_Reference.rc.out"),
        expand(METRICS_DIR / "{run_id}.dedup.metrics.txt",run_id=RUN_IDS),
        expand(ALIGN_DIR / "{run_id}.sorted.filt.bam", run_id=RUN_IDS),
        expand(METRICS_DIR / "{run_id}.filt.metrics.txt", run_id=RUN_IDS), 
        expand(ALIGN_DIR / "{exp_id}.sorted.merged.bam", exp_id=EXPS),
        expand(WIG_DIR / "{exp_id}.coverage.bw", exp_id=EXPS),
        expand(PEAKS_DIR / "{exp_id}_peaks.narrowPeak", exp_id=EXPS),
        expand(PEAKS_DIR / "{exp_id}_peaks.broadPeak", exp_id=EXPS)


# Rule 1. Create BWA Index from FASTA Reference
# -----------------------------------------------------------------------------
bwa_index_rp = cfg.get_rule_params(rulename="BWA_Index_Reference")
rule BWA_Index_Reference:
    input: fasta_path = REF
    params: **(bwa_index_rp.parameters),
        index_dir = REF.parent
    resources: **(bwa_index_rp.resources),
        job_id = "glob",
        logs = LOGS_DIR
    output: rc = (REF.parent / "BWA_Index_Reference.rc.out")
    shell:
        "bwa index {input.fasta_path}"
        " && check_directory -o {output.rc}"
        " {params.index_files} {params.index_dir}"


# Rule 2. Align FASTQs using BWA-MEM algo, then sort, convert to bam, and index 
# -----------------------------------------------------------------------------
bwa_mem_align_rp = cfg.get_rule_params(rulename="BWA_MEM_Align")
rule BWA_MEM_Align:
    input: 
        reads1 = lambda wildcards: get_fastq_file(run_id=f"{wildcards.run_id}", pe_orient=1),
        reads2 = lambda wildcards: get_fastq_file(run_id=f"{wildcards.run_id}", pe_orient=2)
    params: **(bwa_mem_align_rp.parameters),
        ref_prefix = REF
    resources: **(bwa_mem_align_rp.resources),
        job_id = lambda wildcards: f"{wildcards.run_id}",
        logs = LOGS_DIR
    output: 
        bam = temp(SCRATCH_DIR / "{run_id}.sorted.bam"),
        bai = temp(SCRATCH_DIR / "{run_id}.sorted.bam.bai")
    shell:
        "mkdir -p data/align/ &&"
        " bwa mem {params.extra_args} -t {params.threads}"
        " {params.ref_prefix} {input.reads1} {input.reads2}"
        " | samtools sort -@ {params.threads} -o {output.bam} -"
        " && samtools index {output.bam}"


# Rule 3. Mark Duplicates using Picard and remove duplicates
# -----------------------------------------------------------------------------
gatk_markdups_rp = cfg.get_rule_params(rulename="Picard_Mark_Duplicates")
rule Picard_Mark_Duplicates:
    input: 
        bam = rules.BWA_MEM_Align.output.bam,
        bai = rules.BWA_MEM_Align.output.bai
    params: **(gatk_markdups_rp.parameters),
        tmp_dir = lambda wildcards: SCRATCH_DIR / f"{wildcards.run_id}-dedup/temp/"
    resources: **(gatk_markdups_rp.resources), 
        job_id = lambda wildcards: f"{wildcards.run_id}",
        java_mem_max = round(0.75*(gatk_markdups_rp.resources.total_memory_mb)),
        java_mem_min = 2000,
        logs = LOGS_DIR
    output:
        metrics = (METRICS_DIR / "{run_id}.dedup.metrics.txt"),
        bam_dedup = temp(SCRATCH_DIR / "{run_id}-dedup/{run_id}.sorted.dedup.bam"),
        bai_dedup = temp(SCRATCH_DIR / "{run_id}-dedup/{run_id}.sorted.dedup.bam.bai")
    shell:
        "mkdir -p {params.tmp_dir} && "
        "java -Xms{resources.java_mem_min}m -Xmx{resources.java_mem_max}m"
        " -jar $PICARD MarkDuplicates"
        " I={input.bam} M={output.metrics} TMP_DIR={params.tmp_dir}"
        " O={output.bam_dedup} REMOVE_DUPLICATES=true"
        " {params.extra_args} 2> {resources.logs}/{resources.job_id}.MD.log"
        " && samtools index {output.bam_dedup}"


# Rule 4. Shift alignments using deeptools
# -----------------------------------------------------------------------------
deeptools_filt_shift_rp = cfg.get_rule_params(rulename="DeepTools_Filter_And_Shift")
rule DeepTools_Filter_And_Shift:
    input: 
        dedup_bam = rules.Picard_Mark_Duplicates.output.bam_dedup,
        dedup_bai = rules.Picard_Mark_Duplicates.output.bai_dedup,
        blacklist = GLOBALS.files.reference_blacklist
    params: **(deeptools_filt_shift_rp.parameters),
        allowed_chrs = ALLOWED_CHRS
    resources: **(deeptools_filt_shift_rp.resources),
        job_id = lambda wildcards: f"{wildcards.run_id}",
        logs = LOGS_DIR
    output:
        prefilt_bam = temp(SCRATCH_DIR / "{run_id}.sorted.prefilt.bam"),
        prefilt_bai = temp(SCRATCH_DIR / "{run_id}.sorted.prefilt.bam.bai"),
        filt_bam = temp(SCRATCH_DIR / "{run_id}.filt.bam"),
        filt_sort_bam = (ALIGN_DIR / "{run_id}.sorted.filt.bam"),
        filt_sort_bai = (ALIGN_DIR / "{run_id}.sorted.filt.bam.bai"),
        filt_metrics = (METRICS_DIR / "{run_id}.filtshift.metrics.txt")
    shell:
        "samtools view -@ {resources.cpus_per_node}"
        " -b {input.dedup_bam} {params.allowed_chrs}"
        " > {output.prefilt_bam} && samtools index {output.prefilt_bam} && "
        "alignmentSieve -p max -b {output.prefilt_bam} -o {output.filt_bam}"
        " --ATACshift --blackListFileName {input.blacklist}"
        " --filterMetrics {output.filt_metrics} {params.extra_args} && "
        "samtools sort -@ {resources.cpus_per_node} -o {output.filt_sort_bam}"
        " {output.filt_bam} && samtools index {output.filt_sort_bam}"


# Rule 5. Plot diagnostics for individual alignments
# -----------------------------------------------------------------------------
# TODO implement script for fragment-length distribution and TSS-Enrichment

# TODO: implement exclusion lists for excluding problematic samples from merging
# Rule 6. Merge BAMs for each experiment, using Picard and a python wrapper 
# -----------------------------------------------------------------------------
merge_exp_bams_rp = cfg.get_rule_params(rulename="Merge_Experiment_BAMs")
rule Merge_Experiment_BAMs:
    input: bams = expand(ALIGN_DIR / "{run_id}.sorted.filt.bam", run_id=RUN_IDS)
    params: **(merge_exp_bams_rp.parameters), outdir = ALIGN_DIR,
        scripts_dir = SCRIPTS_DIR, run_ids = RUN_IDS,
        exp_ids = EXP_IDS, tmp_dir = (SCRATCH_DIR / "merge_tmp/")
    resources: **(merge_exp_bams_rp.resources),
        java_mem_max = round(0.75*(merge_exp_bams_rp.resources.total_memory_mb)),
        java_mem_min = 2000,
        job_id = "glob"
    output: expand(ALIGN_DIR / "{exp_id}.sorted.merged.bam", exp_id=EXPS),
        expand(ALIGN_DIR / "{exp_id}.bam_merge.list", exp_id=EXPS)
    shell:
        "mkdir -p {params.tmp_dir} && "
        "python {params.scripts_dir}/merge_experiment_bams.py" 
        " --picard-jar $PICARD --bams {input.bams}"
        " --run-ids {params.run_ids} --exp-ids {params.exp_ids}"
        " --outdir {params.outdir}"
        " --picard-extra-args '{params.picard_extra_args} TMP_DIR={params.tmp_dir}'"
        " --java-args '-Xms{resources.java_mem_min}m -Xmx{resources.java_mem_max}m'"


# Rule 7. Rerun diagnostic plots on the merged files
# -----------------------------------------------------------------------------
# TODO: implement scipt for fragment-length distribution and TSS-Enrichment 


# Rule 8. Create bigwig files from Experiment-merged BAMs
# -----------------------------------------------------------------------------
deeptools_bam_coverage_rp = cfg.get_rule_params(rulename="DeepTools_bamCoverage")
rule DeepTools_bamCoverage:
    input: merged_bam = (ALIGN_DIR / "{exp_id}.sorted.merged.bam")
    params: **(deeptools_bam_coverage_rp.parameters),
        wig_dir = WIG_DIR
    resources: **(deeptools_bam_coverage_rp.resources),
        job_id = lambda wildcards: f"{wildcards.exp_id}",
        logs = LOGS_DIR
    output: bw = (WIG_DIR / "{exp_id}.coverage.bw")
    shell:
        "mkdir -p {params.wig_dir} && "
        "bamCoverage --binSize {params.binsize}"
        " --normalizeUsing {params.norm_strategy}"
        " --effectiveGenomeSize {params.effective_gen_size}"
        " --numberOfProcessors {resources.cpus_per_node}"
        " -b {input.merged_bam} -o {output.bw}"


# Rule 9. Call open chromatin peaks using MACS2
# -----------------------------------------------------------------------------
macs2_peakcalling_rp = cfg.get_rule_params(rulename="MACS2_Peak_Calling")
rule MACS2_Peak_Calling:
    input: merged_bam = (ALIGN_DIR / "{exp_id}.sorted.merged.bam")
    params: **(macs2_peakcalling_rp.parameters),
        peaks_dir = PEAKS_DIR
    resources: **(macs2_peakcalling_rp.resources),
        job_id = lambda wildcards: f"{wildcards.exp_id}",
        logs = LOGS_DIR
    output: 
        narrow_peak = (PEAKS_DIR / "{exp_id}_peaks.narrowPeak"),
        broad_peak = (PEAKS_DIR / "{exp_id}_peaks.broadPeak")
    shell:
       "mkdir -p {params.peaks_dir} && "
       "macs2 callpeak -f BAMPE"
       " -t {input.merged_bam} --outdir {params.peaks_dir}"
       " -g hs -B -q {params.qvalue_thresh}"
       " -n {resources.job_id} {params.extra_args} && "
       "macs2 callpeak -f BAMPE"
       " -t {input.merged_bam} --outdir {params.peaks_dir}"
       " -g hs -B -q {params.qvalue_thresh} --broad"
       " -n {resources.job_id} {params.extra_args}"

# Rule 5. Call open chromatin peaks using HOMER
# -----------------------------------------------------------------------------
# TODO implement or delete

# Rule 6. Generate consensus peak sets for cell type
# -----------------------------------------------------------------------------
# TODO implement for experiment




