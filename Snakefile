# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
import os
from snakemake.shell import shell

# --- Defining Some Functions --- #
def checkForGuppyLog(path):
    for dirpath, dirlist, filenames in os.walk(path):
        for name in filenames:
            if name.endswith('.log'):
                return True
    return False

# return run and barcode for given sample ID using
# nested dictonary where values are sampleIDs and
# nested keys are run names
def findSampleFastq(sampleID):
    sample_dict = config['sample_dict']
    for runName in sample_dict:
        for barcode in sample_dict[runName]:
            if sample_dict[runName][barcode] == sampleID:
                return {'runName': runName, 'barcode': barcode}

# --- Importing Configuration File and Defining Important Lists --- #
configfile: "config.yaml"
# Only for pre-basecalled unfiltered guppy reads
MY_RUNNAMES = ["Run0", "Run1_pf_mixed", "Run2_mixed", "Run3_mixed", "Run4_mixed"]
MY_RUNNAMES_QC = ['Exp2_15Nov', 'Exp3_12Dec', 'Exp4_14Mar']


# --- Some rules --- #

rule all:
    input:
        # for basecalling
        # expand(os.path.join("fastq", "{runnames}.fastq"), runnames=config['runnames']),
        # expand(os.path.join("fastq", "{runnames}.fastq"), runnames=MY_RUNNAMES),
        # for runQC
        # expand(os.path.join("QC", "runs", "MinionQC", "{runnamesQC}"), runnamesQC=MY_RUNNAMES_QC),
        # for qcat
        # expand(os.path.join("fastq", "{runnames}.fastq"), runnames=config['runnames']),
        # expand(os.path.join(config['ROOT'], "qcat_trimmed", "{qcat_test_name}"), qcat_test_name=MY_RUNNAMES),
        # accumulate samples
        expand(os.path.join("fastq", "samples", "{samples}.fastq.gz"), samples=config['samples']),
        # for sample QC
        expand(os.path.join("QC", "NanoStat", "{samples}"), samples=config['samples']),
        expand(os.path.join("QC", "NanoPlot", "{samples}"), samples=config['samples'])
    threads: 8


# rule basecalling:
#     input:
#         raw_dir=os.path.join(config['RAWDIR'], "{runnames}/")
#     output:
#         run_fastq=os.path.join("fastq", "{runnames}.fastq")
#     run:
#         # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4) (--min_qscore 7 is already default)
#         # conditionally allow for --resume figure out how later
#         # there is a recent issue April2020 with barcode trimming in guppy_basecaller so use qcat for demult
#         # dna_r9.4.1_450bps_hac.cgf for FLO-MIN106 and SQK-LSK109 combination
#         guppy_output_dir = os.path.join(config['ROOT'], "guppy_output", wildcards.runnames)
#         try:
#             os.makedirs(guppy_output_dir)
#         except FileExistsError:
#             flag = checkForGuppyLog(guppy_output_dir)
#             pass
#         else:
#             flag = False
#         args = {
#         "input":input.raw_dir,
#         "output_dir":guppy_output_dir
#         }
#         command = "guppy_basecaller --resume --input_path {input} --save_path {output_dir} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
#         command = command.format(**args)
#         if flag:
#            shell(command)
#            for dirpath, dirlist, filenames in os.walk(os.path.join(guppy_output_dir, wildcards.runnames)):
#                for name in filenames:
#                    if name.endswith('.fastq'):
#                        os.rename(os.path.join(dirpath, name), os.path.join(dirpath, wildcards.runnames+".fastq"))
#            try:
#                shell("cat "+guppy_output_dir+"/pass/*.fastq > fastq/"+wildcards.runnames+".fastq")
#            except:
#                print("no basecalling happened")
#                shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
#         else:
#             print("No log file to resume from. Starting fresh instance of basecallig")
#             command = "guppy_basecaller --input_path {input} --save_path {output_dir} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
#             command = command.format(**args)
#             shell(command)
#             for dirpath, dirlist, filenames in os.walk(os.path.join(guppy_output_dir, wildcards.runnames)):
#                 for name in filenames:
#                     if name.endswith('.fastq'):
#                         os.rename(os.path.join(dirpath, name), os.path.join(dirpath, wildcards.runnames+".fastq"))
#             try:
#                 shell("cat "+guppy_output_dir+"/pass/*.fastq > fastq/"+wildcards.runnames+".fastq")
#             except:
#                 print("no basecalling happened")
#                 shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
#
# remove empty fastq files to avoid errors later
#
# For run1 no seq summary into output dir copying MinionQC results. Skip Nanocomp QC.
# For run4 considering only first section (> 1/2 of the data) of run other section basecalled as seperate set
# find sequencing summary
# rule runQC:
#     input:
#         seq_summary=os.path.join("guppy_output", "{runnamesQC}", "sequencing_summary.old.txt")
#     output:
#         MinionQC_out=directory(os.path.join("QC", "runs", "MinionQC", "{runnamesQC}")),
#     run:
#         try:
#             os.makedirs(os.path.join("QC", "runs", "MinionQC"))
#         except FileExistsError:
#             pass
#         args = {
#         "input":input.seq_summary,
#         "outputMin":os.path.join("QC", "runs", "MinionQC"),
#         # "minionQCpath":"/media/utlab/DATA_HDD1/Nanopore_metagenomics/Softwares_for_analysis/minion_qc/MinIONQC.R"
#         "minionQCpath":config["minionQCpath"]
#
#         }
#         # shift minionQCpath to config
#         #  -s makes small figures suitable for export rather than optimised for screen
#         command = "Rscript {minionQCpath} -i {input} -o {outputMin} -s TRUE"
#         command = command.format(**args)
#         shell(command)
#
# qcat does trimming simultaneaously if untrimmed files are needed specifically, edit demultiplex_keep_trim
# rule demultiplex_trim:
#     input:
#         raw_fastq="fastq/{qcat_test_name}.fastq"
#     output:
#         trimmed_dir=directory(os.path.join(config['ROOT'], "qcat_trimmed", "{qcat_test_name}"))
#     run:
#         args = {
#         "input":input.raw_fastq,
#         "outputTrimmed":output.trimmed_dir,
#         "kit":config['barcode_kit'],
#         "tsvPath":os.path.join(config['ROOT'], "qcat_trimmed", wildcards.qcat_test_name)
#         }
#         command = "qcat --fastq {input} --barcode_dir {outputTrimmed} --trim -k {kit} --detect-middle --tsv > {tsvPath}.tsv"
#         command = command.format(**args)
#         shell(command)
############################################################################################################
# BELOW RULE DOES NOT WORK AT ALL. IMPLEMENT ONLY IF NEEDED
############################################################################################################
# qcat does trimming simultaneaously so uncomment demultiplex_keep_trim
# if untrimmed files are needed specifically
# BARCODES=[01,02,...]
# rule demultiplex_keep_trim:
#     input:
#         raw_fastq="fastq/{qcat_test_name}.fastq"
#     output:
#         demux_dir=directory(os.path.join(config['ROOT'], "qcat_output/demuxd", "{qcat_test_name}")),
#         trimmed_dir=directory(os.path.join(config['ROOT'], "qcat_output/trimmed", "{qcat_test_name}"))
#     run:
#         args = {
#         "input":os.path.join(output.demux_dir, {qcat_test_name}),
#         # "outputDemux":output.demux_dir,
#         "outputTrimmed":output.trimmed_dir,
#         "kit":config['barcode_kit']
#         }
#         command = "qcat --fastq {input} --barcode_dir {outputTrimmed} -k {kit} --detect-middle"
#         command = command.format(**args)
#         shell(command)
#
#         command = "qcat --fastq {input} --trim -k {kit} --detect-middle"
# Run 0 (sample 01) trimmed in qcat seperately. Will be combined from this point on
############################################################################################################
############################################################################################################
# move to config
# dict of sample name and barcode (called sample_dict) in each run
# X:
#  'barcode': 'sample ID'}
# added to config (X is 0 to 4 for run number)
#
#
rule collectSamples:
    output:
        os.path.join("fastq", "samples", "{samples}.fastq.gz")
    run:
        try:
            os.makedirs(os.path.join("fastq", "samples"))
        except FileExistsError:
            pass
        print(findSampleFastq(wildcards.samples))
        for runName, barcode in findSampleFastq(wildcards.samples).values(), :
            fastqPath = os.path.join("qcat_trimmed", runName, "barcode"+barcode+".fastq")
            if os.path.exists(fastqPath):
                print("\n {} file exists".format(fastqPath))
                shell("cat "+fastqPath+" > {output}")
            else:
                print("\n {} NO file exists".format(fastqPath))
                shell("touch "+fastqPath)
        # zips all files. Unxip as needed for further use
        shell("gzip fastq/samples/*.fastq")
#
# QC of fastq files
rule sampleQC:
    input:
        sampleFastq=os.path.join("fastq", "samples", "{samples}.fastq.gz")
    output:
        nanostat=os.path.join("QC", "NanoStat", "{samples}"),
        nanoplot=directory(os.path.join("QC", "NanoPlot", "{samples}"))
    run:
        args = {
            "input":input.sampleFastq,
            "output_stat":os.path.join("QC", "NanoStat"),
            "output_plot":os.path.join("QC", "NanoPlot", "{samples}"),
            "name":wildcards.samples
        }
        os.makedirs(args['output_plot'])
        shell("touch "+os.path.join("QC", "NanoStat", "{samples}"))
        command_stat = "NanoStat --fastq {input} --outdir {output_stat} -n {name}"
        shell(command_stat.format(**args))
        command_plot = "NanoPlot --verbose --fastq reads.fastq.gz --outdir {output_plot} --prefix {name}"
        shell(command_stat.format(**args))
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
###########--figure out how to do comparative analysis--##########
