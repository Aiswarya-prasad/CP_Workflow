# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
import os

# --- Importing Configuration File --- #
configfile: "config.yaml"


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
rule basecalling:
    input:
        os.path.join(config['RAWDIR'], "{runnames}")
    output:
        directory(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"))
    run:
        os.makedirs(os.path.join(config['RAWDIR'], "guppy_output", wildcards.runnames))
        print(input)
        print(output)

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
# rule basecalling:
#     input:
#         os.path.join(config['RAWDIR'], "{runnames}")
#     output:
#         directory(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"))
#     shell:
#         # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4) (--min_qscore 7 is already default)
#         # conditionally allow for --resume figure out how later
#         # there is a recent issue April2020 with barcode trimming in guppy_basecaller so use qcat for demult
#         # dna_r9.4.1_450bps_hac.cgf for FLO-MIN106 and SQK-LSK109 combination
#         "guppy_basecaller --input_path {input} --save_path {output} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
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
# include run/s that are/were live basecalled or were only available as fastq? USE qcat
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
