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

# find run and barcode for given sample ID using
# returns path to fastq "{runName}/barcode{barcode}.fastq"
def findSampleFastq(wildcards):
    sampleDict = config['sample_dict']
    runBarcodeDict = {}
    for runName in sampleDict:
        for barcode in sampleDict[runName]:
            if sampleDict[runName][barcode] == wildcards:
                runBarcodeDict = {'runName': runName, 'barcode': barcode}
                return join("qcat_trimmed", runBarcodeDict['runName'], "barcode"+runBarcodeDict['barcode']+".fastq")

# can also be a directory with multiple seq summary files for the same input
# modify function accordingly esp for interrupted basecalling
def findSeqSummary(wildcards):
    return join("guppy_output", wildcards, "sequencing_summary.old.txt")

# --- Importing Configuration File and Defining Important Lists --- #
configfile: "config.yaml"
# Only for pre-basecalled unfiltered guppy reads
MY_RUNNAMES = ["Run0", "Run1_pf_mixed", "Run2_mixed", "Run3_mixed", "Run4_mixed"]
MY_RUNNAMES_QC = ['Exp2_15Nov', 'Exp3_12Dec', 'Exp4_14Mar']
GUPPY_RUNNAMES = ['Exp1_25Oct', 'Exp2_15Nov', 'Exp3_12Dec', 'Exp4_14Mar']

# --- Some rules --- #

rule all:
    input:
        #--> for basecalling
        # expand(join("fastq", "{runnames}.fastq"), runnames=config['runnames']),
        expand(join("fastq", "{runnames}.fastq"), runnames=GUPPY_RUNNAMES),
        #--> runQC
        expand(join("QC", "runs", "{runnames}", "{runnames}_NanoStats.txt"), runnames=config['runnames']),
        # expand(join("QC", "runs", "{runnames}", "{runnames}_NanoStats.txt"), runnames=MY_RUNNAMES_QC),
        #--> demultiplex_trim
        expand(join(config['ROOT'], "qcat_trimmed", "{runnames}"), runnames=config['runnames']),
        # expand(join(config['ROOT'], "qcat_trimmed", "{runnames}.tsv"), runnames=config['runnames']),
        # expand(join(config['ROOT'], "qcat_trimmed", "{qcat_test_name}.tsv"), qcat_test_name=MY_RUNNAMES),
        # expand(join(config['ROOT'], "qcat_trimmed", "{qcat_test_name}"), qcat_test_name=MY_RUNNAMES),
        #--> summary
        expand(join(config['ROOT'], "qcat_trimmed", "{runnames}", "summary.txt"), runnames=config['runnames']),
        #--> collectSamples
        # expand(join("fastq", "samples", "{samples}.fastq.gz"), samples=config['samples']),
        #--> sampleQC
        expand(join("QC", "samples", "{samples}", "{samples}_NanoPlot-report.html"), samples=config['samples']),
        #--> kraken2
        expand(join("classified", "{samples}", "kraken2_customdb", "result"), samples=config['samples']),
        expand(join("classified", "{samples}", "kraken2_humandb", "result"), samples=config['samples']),
        expand(join("classified", "{samples}", "kraken2_BacArchViFunProt", "result"), samples=config['samples']),
        #--> bracken
        expand(join("classified", "{samples}", "bracken", "species_report"), samples=config['samples']),
        expand(join("classified", "{samples}", "bracken", "genus_report"), samples=config['samples']),
        #--> centrifuge
        expand(join("classified", "{samples}", "centrifuge", "report"), samples=config['samples'])
    threads: 8


rule basecalling:
    input:
        raw_dir=join(config['RAWDIR'], "{runnames}/")
    output:
        run_fastq=join("fastq", "{runnames}.fastq")
    run:
        # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4) (--min_qscore 7 is already default)
        # conditionally allow for --resume figure out how later
        # there is a recent issue April2020 with barcode trimming in guppy_basecaller so use qcat for demult
        # dna_r9.4.1_450bps_hac.cgf for FLO-MIN106 and SQK-LSK109 combination
        guppy_output_dir = join(config['ROOT'], "guppy_output", wildcards.runnames)
        try:
            makedirs(guppy_output_dir)
        except FileExistsError:
            flag = checkForGuppyLog(guppy_output_dir)
            pass
        else:
            flag = False
        args = {
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
        # try:
        #     makedirs(join("QC", "runs", "MinionQC"))
        # except FileExistsError:
        #     pass
        args = {
        "input":input.seq_summary,
        # "outputMin":join("QC", "runs", "MinionQC"),
        # "minionQCpath":"/media/utlab/DATA_HDD1/Nanopore_metagenomics/Softwares_for_analysis/minion_qc/MinIONQC.R",
        # "minionQCpath":snakemake.config["minionQCpath"]
        # "outputNanoS":join("QC", "runs", "NanoStat"),
        "outputNanoP":join("QC", "runs", wildcards.runnames),
        # "name": wildcards.runnames,
        "prefix": wildcards.runnames+"_"
        }
        # shift minionQCpath to config
        #  -s makes small figures suitable for export rather than optimised for screen
        # command = "Rscript {minionQCpath} -i {input} -o {outputMin} -s TRUE"
        # shell(command.format(**args))
        # command_nanoS = "NanoStat --summary {input} --outdir {outputNanoS} -n {name} --readtype 1D"
        # shell(command_nanoS.format(**args))
        command_nanoP = "NanoPlot --summary {input} --outdir {outputNanoP} -p {prefix} --readtype 1D"
        shell(command_nanoP.format(**args)+" || touch {output}")
#
# qcat does trimming simultaneaously if untrimmed files are needed specifically, edit demultiplex_keep_trim
rule demultiplex_trim:
    input:
        raw_fastq="fastq/{runnames}.fastq"
        # raw_fastq="fastq/{qcat_test_name}.fastq"
    output:
        # trimmed_dir=directory(join(config['ROOT'], "qcat_trimmed", "{qcat_test_name}"))
        trimmed_dir=directory(join(config['ROOT'], "qcat_trimmed", "{runnames}")),
        tsv=join(config['ROOT'], "qcat_trimmed", "{runnames}.tsv")
    run:
        args = {
        "input":input.raw_fastq,
        "outputTrimmed":output.trimmed_dir,
        "kit":config['barcode_kit'],
        # "tsvPath":join(config['ROOT'], "qcat_trimmed", wildcards.qcat_test_name)
        "tsvPath":join(config['ROOT'], "qcat_trimmed", wildcards.runnames)
        }
        command = "qcat --fastq {input} --barcode_dir {outputTrimmed} --trim -k {kit} --detect-middle --tsv > {tsvPath}.tsv"
        command = command.format(**args)
        shell(command)

rule demultiplex_summary:
    input:
        tsv=join(config['ROOT'], "qcat_trimmed", "{runnames}.tsv")
    output:
        txt=join(config['ROOT'], "qcat_trimmed", "{runnames}", "summary.txt"),
        png=join(config['ROOT'], "qcat_trimmed", "{runnames}", "summary.png")
    script:
        "scripts/demultiplex_summarize.py"


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
#         demux_dir=directory(join(config['ROOT'], "qcat_output/demuxd", "{qcat_test_name}")),
#         trimmed_dir=directory(join(config['ROOT'], "qcat_output/trimmed", "{qcat_test_name}"))
#     run:
#         args = {
#         "input":join(output.demux_dir, {qcat_test_name}),
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
# dict of sample name and barcode (called sampleDict) in each run
# X:
#  'barcode': 'sample ID'}
# added to config (X is 0 to 4 for run number)
#
#
rule collectSamples:
    input:
        # fastqPath=join(config['ROOT'], "qcat_trimmed", findSampleFastq("{samples}"))
        fastqPath=lambda wildcards: findSampleFastq(wildcards.samples)
    output:
        # join("fastq", "samples", "{runnames}_{samples}.fastq.gz")
        join("fastq", "samples", "{samples}.fastq.gz")
    run:
        try:
            makedirs(join("fastq", "samples"))
        except FileExistsError:
            pass
        if exists(input.fastqPath):
            print("\n {} file exists".format(input.fastqPath))
            shell("cat "+input.fastqPath+" > "+join("fastq", "samples", wildcards.samples+".fastq"))
        else:
            print("\n {} NO file exists".format(input.fastqPath))
            shell("touch "+input.fastqPath)
        # zips all files. Unxip as needed for further use
        shell("gzip "+join("fastq", "samples", wildcards.samples+".fastq"))
#
# QC of fastq files
rule sampleQC:
    input:
        sampleFastq=join("fastq", "samples", "{samples}.fastq.gz")
    output:
        # nanostat=join("QC", "NanoStat", "{samples}"),
        # nanoplot=directory(join("QC", "NanoPlot", "{samples}"))
        Nanoplot_Dynamic_Histogram_Read_length_html = join("QC", "samples", "{samples}", "{samples}_Dynamic_Histogram_Read_length.html"),
        Nanoplot_HistogramReadlength_png = join("QC", "samples", "{samples}", "{samples}_HistogramReadlength.png"),
        Nanoplot_LengthvsQualityScatterPlot_dot_png = join("QC", "samples", "{samples}", "{samples}_LengthvsQualityScatterPlot_dot.png"),
        Nanoplot_LengthvsQualityScatterPlot_kde_png = join("QC", "samples", "{samples}", "{samples}_LengthvsQualityScatterPlot_kde.png"),
        Nanoplot_LogTransformed_HistogramReadlength_png = join("QC", "samples", "{samples}", "{samples}_LogTransformed_HistogramReadlength.png"),
        Nanoplot_report_html = join("QC", "samples", "{samples}", "{samples}_NanoPlot-report.html"),
        Nanoplot_NanoStats_txt = join("QC", "samples", "{samples}", "{samples}_NanoStats.txt"),
        Nanoplot_Weighted_HistogramReadlength_png = join("QC", "samples", "{samples}", "{samples}_Weighted_HistogramReadlength.png"),
        Nanoplot_Weighted_LogTransformed_HistogramReadlength_png = join("QC", "samples", "{samples}", "{samples}_Weighted_LogTransformed_HistogramReadlength.png"),
        Nanoplot_Yield_By_Length_png = join("QC", "samples", "{samples}", "{samples}_Yield_By_Length.png")
    run:
        args = {
            "input":input.sampleFastq,
            # "output_stat":join("QC", "NanoStat"),
            "output_plot":join("QC", "samples", wildcards.samples),
            "name":wildcards.samples,
            "prefix":wildcards.samples+'_'
        }
        command = "NanoPlot --verbose --fastq {input} --outdir {output_plot} --prefix {prefix}"
        shell(command.format(**args))
#
# include run/s that are/were live basecalled or were only available as fastq?
#
# "C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
# The sequence ID, obtained from the FASTA/FASTQ header.
# The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.
# The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".
# A space-delimited list indicating the LCA mapping of each k-mer in the sequence(s). For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
# the first 13 k-mers mapped to taxonomy ID #562
# the next 4 k-mers mapped to taxonomy ID #561
# the next 31 k-mers contained an ambiguous nucleotide
# the next k-mer was not in the database
# the last 3 k-mers mapped to taxonomy ID #562
# rule kraken2_human:
#     input:
#         fastq=join("fastq", "samples", "{samples}.fastq.gz")
#     output:
#         report_mpa=join("classified", "{samples}", "kraken2_customdb", "report_mpa"),
#         report_kraken=join("classified", "{samples}", "kraken2_customdb", "report"),
#         result=join("classified", "{samples}", "kraken2_customdb", "result"),
#         unclass=join("classified", "{samples}", "humanDB_unclassified")
#     run:
#         args = {
#         "db": config['kraken_db'],
#         "t": 8,
#         "input": input.fastq,
#         "output_reportmpa": output.report_mpa,
#         "output_report": output.report_kraken,
#         "output_result": output.result,
#         "unclass_out": output.unclass
#         }
#         commandMPA = "kraken2 --db {db} --threads {t}  --gzip-compressed {input} --report {output_reportmpa} --report-zero-counts --use-mpa-style --output {output_result}"
#         shell(commandMPA.format(**args))
#         command = "kraken2 --db {db} --threads {t}  --gzip-compressed {input} --report {output_report} --report-zero-counts --output {output_result} --unclassified-out {unclass_out}"
#         shell(command.format(**args))

rule kraken2:
    input:
        fastq=join("fastq", "samples", "{samples}.fastq.gz")
    output:
        report_mpa=join("classified", "{samples}", "kraken2_customdb", "report_mpa"),
        report_krakendb=join("classified", "{samples}", "kraken2_customdb", "report"),
        result_krakendb=join("classified", "{samples}", "kraken2_customdb", "result"),
        report_humandb=join("classified", "{samples}", "kraken2_humandb", "report"),
        result_humandb=join("classified", "{samples}", "kraken2_humandb", "result"),
        report_customdb=join("classified", "{samples}", "kraken2_BacArchViFunProt`", "report"),
        result_customdb=join("classified", "{samples}", "kraken2_BacArchViFunProt", "result")
    run:
        args = {
        "db": config['kraken_db'],
        "db_human": join("media", "utlab", "DATA_HDD1", "Nanopore_metagenomics", "Softwares_for_analysis", "kraken2", "dbs", "db_humanVec_May2020"),
        "db_Bac": join("media", "utlab", "DATA_HDD1", "Nanopore_metagenomics", "Softwares_for_analysis", "kraken2", "dbs", "db_BacArchViFunProt_May2020"),
        "t": 8,
        "input": input.fastq,
        "output_reportmpa": output.report_mpa,
        "output_report_krakendb": output.report_krakendb,
        "output_result_krakendb": output.result_krakendb,
        "output_report_human": output.report_human,
        "output_result_human": output.result_human,
        "output_report_custom": output.report_custom,
        "output_result_custom": output.result_custom
        }
        commandMPA = "kraken2 --db {db} --threads {t}  --gzip-compressed {input} --report {output_reportmpa} --report-zero-counts --use-mpa-style --output {output_result}"
        command = "kraken2 --db {db} --threads {t}  --gzip-compressed {input} --report {output_report_krakendb} --report-zero-counts --output {output_result_krakendb}"
        command_human = "kraken2 --db {db_human} --threads {t}  --gzip-compressed {input} --report {output_report_human} --report-zero-counts --output {output_result_human}"
        command_custom = "kraken2 --db {db_custom} --threads {t}  --gzip-compressed {input} --report {output_report_custom} --report-zero-counts --output {output_result_custom}"
        shell(commandMPA.format(**args))
        shell(command.format(**args))
#
#
rule bracken:
    input:
        kraken_report=rules.kraken2.output.report_kraken
    output:
        reportS=join("classified", "{samples}", "bracken_customdb", "species_report"),
        reportG=join("classified", "{samples}", "bracken_customdb", "genus_report")
    run:
        args = {
        "db": config['kraken_db'],
        "input": input.kraken_report,
        "output_reportS": output.reportS,
        "output_reportG": output.reportG
        }
        commandS = "bracken -d {db} -i {input} -l S -o {output_reportS}"
        commandG = "bracken -d {db} -i {input} -l G -o {output_reportG}"
        shell(commandS.format(**args))
        shell(commandG.format(**args))
#
#
rule centrifuge:
    input:
        fastq=join("fastq", "samples", "{samples}.fastq.gz")
    output:
        report=join("classified", "{samples}", "centrifuge", "report"),
        result=join("classified", "{samples}", "centrifuge", "result")
    run:
        args = {
        "input": input.fastq,
        "db": config['centrifuge_db'],
        "output_report": output.report,
        "output_result": output.result
        }
        # -U - means take from stdin
        command = "gunzip -c {input} | centrifuge -x {db} -q  -U - --report-file {output_report} -S {output_result}"
        shell(command.format(**args))
#
###########--comparative analysis--##########
# rule compare:
#     input:
#         reports, results
#     output:
#         figure
#     script:
#         r
#         py
