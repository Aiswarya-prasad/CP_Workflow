# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
import os

# --- Importing Configuration File --- #
configfile: "config.yaml"

def FindFast5(run_path_list):
    """
    when provided with path to directory containing the
    raw data, it returns the path to the fast5 folder
    """
    fast5_path_list=[]
    for run_path in run_path_list:
        for dirpath, dirnames, filenames in os.walk(run_path):
            if os.path.basename(dirpath) == "fast5":
                for file in os.listdir(dirpath):
                    if os.path.splitext(file)[1] == ".fast5":
                        fast5_path_list.append(dirpath)
                    else:
                        print(dirpath, "contains NO fast5 files")
    return fast5_path_list
# RAWDIR="/media/utlab/DATA"
RAWDIR="/home/aiswarya/UtpalTatu_MS/Metagenomics/CP_Workflow/test"
# RUNNAMES=['Exp1_25Oct', 'Exp2_15Nov', 'Exp3_12Dec', 'Exp4_14Mar']
RUNNAMES=['temp1', 'temp2', 'temp3', 'temp4']
SAMPLES=['01', '03', '06', '07', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22']

# --- Some rules --- #

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
rule basecalling:
    input:
        FindFast5(expand(os.path.join(RAWDIR, "{run_names}"), run_names=RUNNAMES))
    # output:
    #     "consider syntax and documentation"
    shell:
        "echo {input}"
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
