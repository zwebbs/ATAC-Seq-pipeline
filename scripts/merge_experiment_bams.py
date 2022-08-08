# File Name: merge_experiment_bams.py
# Created By: ZW
# Created On: 2022-08-07
# Purpose: processs experiment metadata and merge bam files from each
# experiment into a single alignment set. This is designed for technical
# replicate merging.


# Module + Library Imports
# -----------------------------------------------------------------------------
import argparse
from pathlib import Path
import subprocess
from sys import stderr

# Function Definitions
# -----------------------------------------------------------------------------

# Main execution
#------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # define and implement cmdline parser for passing args to script
    descr = """
     File Name: merge_experiment_bams.py
     Created By: ZW
     Created On: 2022-08-07
     Purpose: processs experiment metadata and merge bam files from each
      experiment into a single alignment set. This is designed for technical
      replicate merging.
    """
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("--bams", type=str,nargs='+',
            help="individual bam files to merge")
    parser.add_argument("--run-ids", type=str, nargs='+',
            help="run ids corresponding to bam files")
    parser.add_argument("--exp-ids",type=str, nargs='+',
            help="string representing sequence map dictionary")
    parser.add_argument("--outdir", type=str,
            help="directory to output merge.list files")
    parser.add_argument("--picard-jar", type=str,
            help="path to the picard jar file")
    parser.add_argument("--java-args", type=str, default="",
            help="string representing the java runtime arguments")
    parser.add_argument("--picard-extra-args", type=str, default="",
            help="string representing extra Picard arguments")
    args = parser.parse_args()
    
    # validate output directory
    outdir = Path(args.outdir).resolve()
    if not outdir.exists():
        parser.error("Error! Path given to --outdir must exist")
        exit(1)

    # get unique experiment ids and output logs
    exps = list(set(args.exp_ids))
    merge_list_filepaths = [outdir / f"{exp}.bam_merge.list" for exp in exps]
    exp_info = zip(exps, merge_list_filepaths)
    
    # get run level information

    # Execute merge routine using subprocess and Picard
    for experiment_name, output_list in exp_info:
        run_ids, exp_ids, bams = [],[],[]
        run_info = zip(args.run_ids, args.exp_ids, args.bams)
        
        # iterate through run info
        for r, e, b in run_info:
            if e == experiment_name:
                run_ids.append(r)
                exp_ids.append(e)
                bams.append(b)
        
        # do not proceed if there are no outputs for the given experiment
        if any([len(run_ids) == 0, len(exp_ids) == 0, len(bams) == 0]):
            continue

        # build Picard command
        input_string = 'I=' + ' I='.join(bams)
        output_string = 'O=' + str(outdir / f"{experiment_name}.sorted.merged.bam")
        
        cmd = (f"java {args.java_args} -jar {args.picard_jar} MergeSamFiles" + 
                f" {input_string} {output_string} {args.picard_extra_args}")
        
        # execute comand and print output logs
        ps = subprocess.run(cmd, capture_output=True, shell=True)

        print("------------------------------")
        print(f"EXP: {experiment_name} -STDOUT")
        print("------------------------------")
        print(f"Running Command: {cmd}\n") 
        print(ps.stdout.decode('utf-8'))
        print("------------------------------", file=stderr)
        print(f"EXP: {experiment_name} -STDERR", file=stderr)
        print("------------------------------", file=stderr)
        print(f"Running Command: {cmd}\n")
        print(ps.stderr.decode('utf-8'), file=stderr)

        # write the log
        with open(output_list, 'w') as outobj:
            run_info = zip(run_ids,exp_ids,bams)
            for r,e,b in run_info:
                outobj.write(f'{r}\t{e}\t{b}\n')



    
