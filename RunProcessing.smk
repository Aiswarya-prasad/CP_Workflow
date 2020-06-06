# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
from os import walk, rename, makedirs
from os.path import join, exists
from snakemake.shell import shell
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import pandas as pd

# --- Defining Some Functions --- #
def checkForGuppyLog(path):
    for dirpath, dirlist, filenames in walk(path):
        for name in filenames:
            if name.endswith('.log'):
                return True
    return False
#
# can also be a directory with multiple seq summary files for the same input
# modify function accordingly esp for interrupted basecalling
def findSeqSummary(wildcards):
    return join("guppy_output", wildcards, "sequencing_summary.old.txt")

# --- Importing Configuration File and Defining Important Lists --- #
configfile: "config.yaml"

# --- The rules --- #

# EDIT THIS AS NEEDED
wildcard_constraints:
    runnames = "^Run*",
    samples = "\d+"

rule all:
    input:
        #--> for basecalling
        # expand(join("fastq", "{runnames}.fastq"), runnames=config['runnames']),
        #--> runQC
        expand(join("QC", "runs", "{runnames}", "{runnames}_NanoStats.txt"), runnames=config['runnames']),
        #--> demultiplex_trim
        expand(join("qcat_trimmed", "{runnames}.tsv"), runnames=config['runnames']),
        expand(join("qcat_trimmed", "{runnames}"), runnames=config['runnames']),
        #--> summary
        expand(join("qcat_trimmed", "{runnames}", "summary.txt"), runnames=config['runnames']),
        expand(join("qcat_trimmed", "{runnames}", "summary.png"), runnames=config['runnames']),
        #--> collectSamples (Do not have to specifi because other rules depend on this)


# edit input as needed by guppy. Current method works for input which is a symlink to the
# since --recursive is used it will search subtree of input which needs to be a directory (or link)
# directory containing the fast5 directory in its subtree
# Also edit config['RAWDIR as needed']
rule basecalling:
    input:
        # raw_dir=join(config['RAWDIR'], "{runnames}") # remove extra / (depending on Rawdir structure)
        raw_dir=join("RawDir", "{runnames}")
    output:
        runFastq=join("fastq", "{runnames}.fastq")
    run:
        # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4) (--min_qscore 7 is already default)
        # conditionally allow for --resume figure out how later
        # there is a recent issue April2020 with barcode trimming in guppy_basecaller so use qcat for demult
        # dna_r9.4.1_450bps_hac.cgf for FLO-MIN106 and SQK-LSK109 combination
        guppy_output_dir = join("guppy_output", wildcards.runnames)
        try:
            makedirs(guppy_output_dir)
        except FileExistsError:
            flag = checkForGuppyLog(guppy_output_dir)
            pass
        else:
            flag = False
        args = {
        # path to dir containing input w/ wildcard
        "input":input.raw_dir,
        "output_dir":guppy_output_dir
        }
        command = "guppy_basecaller --resume --input_path {input} --save_path {output_dir} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
        command = command.format(**args)
        if flag:
           shell(command)
           for dirpath, dirlist, filenames in walk(join(guppy_output_dir, wildcards.runnames)):
               for name in filenames:
                   if name.endswith('.fastq'):
                       rename(join(dirpath, name), join(dirpath, wildcards.runnames+".fastq"))
           try:
               shell("cat "+guppy_output_dir+"/pass/*.fastq > fastq/"+wildcards.runnames+".fastq")
           except:
               print("no basecalling happened")
               shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
        else:
            print("No log file to resume from. Starting fresh instance of basecallig")
            command = "guppy_basecaller --input_path {input} --save_path {output_dir} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
            command = command.format(**args)
            shell(command)
            for dirpath, dirlist, filenames in walk(join(guppy_output_dir, wildcards.runnames)):
                for name in filenames:
                    if name.endswith('.fastq'):
                        rename(join(dirpath, name), join(dirpath, wildcards.runnames+".fastq"))
            try:
                shell("cat "+guppy_output_dir+"/pass/*.fastq > fastq/"+wildcards.runnames+".fastq")
            except:
                print("no basecalling happened")
                shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
#
# remove empty fastq files to avoid errors later
#
# For run1 no seq summary into output dir copying MinionQC results. Skip Nanocomp QC.
# For run4 considering only first section (> 1/2 of the data) of run other section basecalled as seperate set
# find sequencing summary
rule runQC:
    input:
        seq_summary=lambda wildcards: findSeqSummary(wildcards.runnames)
    output:
        # MinionQC_out=directory(join("QC", "runs", "MinionQC", "{runnames}")),
        # NanoStat_out=join("QC", "runs", "NanoStat", "{runnames}"),
        Nanoplot_Dynamic_Histogram_Read_length_html = join("QC", "runs", "{runnames}", "{runnames}_Dynamic_Histogram_Read_length.html"),
        Nanoplot_HistogramReadlength_png = join("QC", "runs", "{runnames}", "{runnames}_HistogramReadlength.png"),
        Nanoplot_LengthvsQualityScatterPlot_dot_png = join("QC", "runs", "{runnames}", "{runnames}_LengthvsQualityScatterPlot_dot.png"),
        Nanoplot_LengthvsQualityScatterPlot_kde_png = join("QC", "runs", "{runnames}", "{runnames}_LengthvsQualityScatterPlot_kde.png"),
        Nanoplot_LogTransformed_HistogramReadlength_png = join("QC", "runs", "{runnames}", "{runnames}_LogTransformed_HistogramReadlength.png"),
        Nanoplot_report_html = join("QC", "runs", "{runnames}", "{runnames}_NanoPlot-report.html"),
        Nanoplot_NanoStats_txt = join("QC", "runs", "{runnames}", "{runnames}_NanoStats.txt"),
        Nanoplot_Weighted_HistogramReadlength_png = join("QC", "runs", "{runnames}", "{runnames}_Weighted_HistogramReadlength.png"),
        Nanoplot_Weighted_LogTransformed_HistogramReadlength_png = join("QC", "runs", "{runnames}", "{runnames}_Weighted_LogTransformed_HistogramReadlength.png"),
        Nanoplot_Yield_By_Length_png = join("QC", "runs", "{runnames}", "{runnames}_Yield_By_Length.png")
    run:
        args = {
        "input":input.seq_summary,
        "outputNanoP":join("QC", "runs", wildcards.runnames),
        "prefix": wildcards.runnames+"_"
        }
        command_nanoP = "NanoPlot --summary {input} --outdir {outputNanoP} -p {prefix} --readtype 1D"
        shell(command_nanoP.format(**args)+" || touch {output}")
#
# qcat does trimming simultaneaously if untrimmed files are needed specifically, edit demultiplex_keep_trim
rule demultiplexTrim:
    input:
        raw_fastq=rules.basecalling.output.runFastq
    output:
        outDir=directory(join("qcat_trimmed", "{runnames}")),
        tsv=join("qcat_trimmed", "{runnames}.tsv")
    run:
        args = {
            "input":input.raw_fastq,
            "outputTrimmed":join(config['ROOT'], "qcat_trimmed", wildcards.runnames),
            "kit":config['barcode_kit'],
            "tsvPath":join(config['ROOT'], "qcat_trimmed", wildcards.runnames)
        }
        command = "qcat --fastq {input} --barcode_dir {outputTrimmed} --trim -k {kit} --detect-middle --tsv > {tsvPath}.tsv"
        command = command.format(**args)
        shell(command)

rule demultiplexSummary:
    input:
        demuxDirs=rules.demultiplexTrim.output.tsv
    output:
        txt=join("qcat_trimmed", "{runnames}", "summary.txt"),
        png=join("qcat_trimmed", "{runnames}", "summary.png")
    script:
        "scripts/demultiplex_summarize.py"
#
#
rule collectSamples:
    input:
        demuxDirs=rules.demultiplexTrim.output.outDir
    output:
        fastq=expand(join("fastq", "samples", "{samples}.fastq.gz"), samples=config['samples']),
    run:
        try:
            makedirs(join("fastq", "samples"))
        except FileExistsError:
            pass
        sampleDict = config['sample_dict']
        runBarcodeDict = {}
        for runName in wildcards.runnames:
            for barcode in sampleDict[runName]:
                if sampleDict[runName][barcode] in config['sample']:
                    sampleName = sampleDict[runName][barcode]
                    runBarcodeDict = {'runName': runName, 'barcode': barcode}
                    return join("qcat_trimmed", runBarcodeDict['runName'], "barcode"+runBarcodeDict['barcode']+".fastq")
                    fastqPath = join("qcat_trimmed", runBarcodeDict['runName'], "barcode"+runBarcodeDict['barcode']+".fastq")
                    command = 'cat '+fastqPath+' | gzip > '+join("fastq", "samples", sampleName+".fastq.gz")
