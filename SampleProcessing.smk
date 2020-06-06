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


# --- Importing Configuration File and Defining Important Lists --- #
configfile: "config.yaml"

# --- The rules --- #

# EDIT THIS AS NEEDED
wildcard_constraints:
    runnames = "^Run*",
    samples = "\d+"

subworkflow RunProcessing:
    snakefile:
        "RunProcessing.smk"

rule all:
    input:
        #--> collectSamples (Do not have to specifi because other rules depend on this)
        RunProcessing(expand(join("fastq", "samples", "{samples}.fastq.gz"), samples=config['samples'])),
        #--> sampleQC
        expand(join("QC", "samples", "{samples}", "{samples}_NanoPlot-report.html"), samples=config['samples']),
        # #--> kraken2
        expand(join("classified", "{samples}", "kraken2_Minidb", "result"), samples=config['samples']),
        # # uncomment if using
        # # expand(join("classified", "{samples}", "kraken2_humandb", "result"), samples=config['samples']),
        # # expand(join("classified", "{samples}", "kraken2_custom", "result"), samples=config['samples']),
        # #--> bracken (depends on Kraken2 report)
        expand(join("classified", "{samples}", "bracken", "species_report"), samples=config['samples']),
        expand(join("classified", "{samples}", "bracken", "genus_report"), samples=config['samples']),
        # #--> centrifuge
        expand(join("classified", "{samples}", "centrifuge", "report"), samples=config['samples']),
        expand(join("classified", "{samples}", "centrifuge", "result"), samples=config['samples'])
#
# QC of fastq files
rule sampleQC:
    input:
        sampleFastq=join("fastq", "samples", "{samples}.fastq.gz"),
    output:
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
            "output_plot":join("QC", "samples", wildcards.samples),
            "name":wildcards.samples,
            "prefix":wildcards.samples+'_'
        }
        command = "NanoPlot --verbose --fastq {input} --outdir {output_plot} --prefix {prefix}"
        shell(command.format(**args))
#
#
rule filterSamples:
    input:
        input=join("fastq", "samples", "{samples}.fastq.gz")
    output:
        output7=join("fastq", "samplesQ7", "{samples}.fastq.gz"),
        output10=join("fastq", "samplesQ10", "{samples}.fastq.gz")
    run:
        args = {
        "output7":output.output7,
        "output10":output.output10
        }
        makedirs(output.output7)
        makedirs(output.output10)
        command7 = "gunzip -c {input} | NanoFilt --quality 7 | gzip > {output7}"
        shell(command7.format(**args))
        command10 = "gunzip -c {input} | NanoFilt --quality 10 | gzip > {output7}"
        shell(command10.format(**args))
#
#
rule kraken2:
    input:
        fastq=join("fastq", "samplesQ7", "{samples}.fastq.gz")
    output:
        # report_mpa=join("classified", "{samples}", "kraken2_Minidb", "report_mpa"),
        report=join("classified", "{samples}", "kraken2_Minidb", "report"),
        result=join("classified", "{samples}", "kraken2_Minidb", "result"),
    run:
        args = {
        "db": config['kraken_db'],
        "t": 8,
        "input": input.fastq,
        # "output_reportmpa": output.report_mpa,
        "output_report": output.report,
        "output_result": output.result,
        }
        commandMPA = "kraken2 --db {db} --confidence ? --threads {t}  --gzip-compressed {input} --report {output_reportmpa} --report-zero-counts --use-mpa-style > /dev/null"
        command = "kraken2 --db {db} --confidence 0.00001 --threads {t}  --gzip-compressed {input} --report {output_report} --report-zero-counts --output {output_result}"
        # shell(commandMPA.format(**args))
        shell(command.format(**args))
#
#
rule bracken:
    input:
        kraken_report=rules.kraken2.output.report
    output:
        reportS=join("classified", "{samples}", "bracken", "species_report"),
        reportG=join("classified", "{samples}", "bracken", "genus_report")
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
        fastq=join("fastq", "samplesQ7", "{samples}.fastq.gz")
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
#
