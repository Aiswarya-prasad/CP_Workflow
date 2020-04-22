# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
import os

# --- Importing Configuration File --- #
configfile: "config.yaml"

# def FindFast5(run_path_list):
#     """
#     when provided with a list of directory paths with
#     raw data (per run), it returns a list of paths to
#     the fast5 directories in them which contain fast5 files
#
#     it can handle multiple fast5 directories in each
#     parent given
#     the path will contain the parent directory (run name)
#     """
#     # for root, dirs, files in os.walk(run_path, topdown=False):
#     #     for name in dirs:
#     #         if os.path.basename(dirpath) == "fast5":
#     #             print(name)
#
#     fast5_path_list=[]
#     for run_path in run_path_list:
#         for dirpath, dirnames, filenames in os.walk(run_path, topdown = False):
#             if os.path.basename(dirpath) == "fast5":
#                 for file in os.listdir(dirpath):
#                     if os.path.splitext(file)[1] == ".fast5":
#                         fast5_path_list.append(dirpath)
#                     else:
#                         print(dirpath, "contains NO fast5 files")
#     return fast5_path_list



# --- Some rules --- #

rule all:
    input:
        expand(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"), runnames=config['runnames'])
threads:2

# JUAT 2 RULES THAT WORK
# rule all:
#     input:
#         expand(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"), runnames=RUNNAMES)
# threads:2
#
# rule basecalling:
#     input:
#         os.path.join(config['RAWDIR'], "{runnames}")
#     output:
#         directory(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"))
#     run:
#         os.makedirs(os.path.join(config['RAWDIR'], "guppy_output", wildcards.runnames))
#         print(input)
#         print(output)

# rule parse_metadata:
    # somehow parse metadata and make it available to all the other rules as needed
    # do this in config?

#
# rule raw_QC:
#     input:
#         "path to sequencing_summary.txt"
#     output:
#         "figures",
#         "text"
#     shell:
#         "program for doing that minionQC (?????)"
#
#         ######## --- will depend on --- ########
#         rule find_sumary:
#             "path to Run dir"
#         output:
#             "path to seq summary"
#         script:
#             "scripts/find_sumary.py"
#
# make sure to consider multiple fast5 files per run in following steps
rule basecalling:
    input:
        os.path.join(config['RAWDIR'], "{runnames}")
    output:
        directory(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"))
    shell:
        # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4)
        "guppy_basecaller --input_path {input} --save_path {output} --config guppy_config.cfg --recursive --records_per_fastq 0 --calib"
#
# rule fastq_QC_run:
#     input:
#         "basecalled fastq files per run"
#     output:
#         "figures",
#         "text or csv"
#     shell:
#         "appropriate program to do QC. If minionQC works, drop this"
#
# include run/s that are/were live basecalled or were only available as fastq?
#
# rule demultiplex:
#     input:
#         "basecalled fastq files per run"
#     output:
#         "all barcodes for each run"
#     shell:
#         "guppy_barcoder (ok? -- other options)"
#
#         ######## --- will depend on --- ########
#         # sample, run, barcode correlation from config or metadata (latter better)
#         rule seggregate_fastq:
#             input:
#                 "result of demux per run"
#             output:
#                 "fastq for each sample"
#             script:
#                 "script/seggregate_fastq.py"
#
# include run/s that are/were live basecalled or were only available as fastq?
#
# rule run_kraken2:
#     input:
#         "fatq files per sample"
#     output:
#         "reports",
#         "results"
#     shell:
#         "call run_kraken2"
#
# rule run_centrifuge:
#     input:
#         "fatq files per sample"
#     output:
#         "reports",
#         "results"
#     shell:
#         "call run_centrifuge"
#
# ###########--figure out how to do comparative analysis--###########
