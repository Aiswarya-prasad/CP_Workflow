# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Some Packages --- #
from os import walk, rename, makedirs
from os.path import join, exists, getmtime
from snakemake.shell import shell
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import pandas as pd

# --- Defining Some Functions --- #


# --- Importing Configuration File and Defining Important Lists --- #
configfile: "config.yaml"

# depends on type of barcode kit used
# BARCODES = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

# --- The rules --- #

rule complete:
    input:
        expand(join("QC", "samples", "{samples}", "{samples}_NanoPlot-report.html"), samples=config['samples']),
        expand(join("00_RawData", "samples_Q7", "{samples}.fastq.gz"), samples=config['samples']),
        expand(join("00_RawData", "samples_Q10", "{samples}.fastq.gz"), samples=config['samples'])
    threads: 8

rule sampleQC:
    input:
        sampleFastq=join("00_RawData", "samples", "{samples}.fastq.gz")
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
        fastq=join("00_RawData", "samples", "{samples}.fastq.gz")
    output:
        output7=join("00_RawData", "samples_Q7", "{samples}.fastq.gz"),
        output10=join("00_RawData", "samples_Q10", "{samples}.fastq.gz")
    run:
        args = {
        "input": input.fastq,
        "output7":output.output7,
        "output10":output.output10
        }
        try:
            makedirs(join("00_RawData", "samples_Q7"))
            makedirs(join("00_RawData", "samples_Q7"))
            makedirs(join("00_RawData", "samples_Q10"))
            makedirs(join("00_RawData", "samples_Q10"))
        except:
            pass
        command7 = "gunzip -c {input} | NanoFilt --quality 7 | gzip > {output7}"
        shell(command7.format(**args))
        command10 = "gunzip -c {input} | NanoFilt --quality 10 | gzip > {output10}"
        shell(command10.format(**args))
#
#
rule unzip:
  input:
    gz=join("00_RawData", "samples", "{samples}.fastq.gz"),
    gz7=join("00_RawData", "samples_Q7", "{samples}.fastq.gz"),
    gz10=join("00_RawData", "samples_Q10", "{samples}.fastq.gz")
  output:
    fq=join("00_RawData", "samples", "{samples}.fastq"),
    fq7=join("00_RawData", "samples_Q7", "{samples}.fastq"),
    fq10=join("00_RawData", "samples_Q10", "{samples}.fastq")
  run:
    shell("gunzip -c "+input.gz+" > "+output.fq)
    # uncoment below line if you want to remove the zipped
    # shell("rm -rf "+input.gz)
    shell("gunzip -c "+input.gz7+" > "+output.fq7)
    # uncoment below line if you want to remove the zipped
    # shell("rm -rf "+input.gz)
    shell("gunzip -c "+input.gz10+" > "+output.fq10)
    # uncoment below line if you want to remove the zipped
    # shell("rm -rf "+input.gz)
#
#
# rule assemble:
#     input:
#         fastq=("fastq", "{runnames}.fastq")
#     output:
#         findout
#     run:
#         flye
