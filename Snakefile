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


# --- Importing Configuration File --- #
configfile: "config.yaml"
# Only for pre-basecalled unfiltered guppy reads
MY_RUNNAMES = ["Run0", "Run1_pf_mixed", "Run2_mixed", "Run3_mixed", "Run4_mixed"]
MY_RUNNAMES_QC = ['Exp2_15Nov', 'Exp3_12Dec', 'Exp4_14Mar']

# --- Some rules --- #

rule all:
    input:
        # expand(os.path.join("fastq", "{runnames}.fastq"), runnames=config['runnames']),
        # expand(os.path.join("QC", "runs", "MinionQC", "{runnames}"), runnames=config['runnames'])
        # expand(os.path.join("QC", "runs", "MinionQC"))
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
#                shell("rsync -v "+guppy_output_dir+"/pass/"+wildcards.runnames+".fastq fastq/"+wildcards.runnames+".fastq")
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
#                 shell("rsync -v "+guppy_output_dir+"/pass/"+wildcards.runnames+".fastq fastq/"+wildcards.runnames+".fastq")
#             except:
#                 print("no basecalling happened")
#                 shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
#
# remove empty fastq files to avoid errors later
#
# For run1 no seq summary into output dir copying MinionQC results. Skip Nanocomp QC.
# For run4 considering only first section (> 1/2 of the data) of run other section basecalled as seperate set
# find sequencing summary
rule runQC:
    input:
        # can also be a directory with multiple for the same input. Use if basecalling was resumed.
        seq_summary=os.path.join("guppy_output", "{runnames}", "sequencing_summary.old.txt")
    output:
        # MinionQC_out=directory(os.path.join("QC", "runs", "MinionQC", "{runnames}")),
        NanoStat_out=os.path.join("QC", "runs", "NanoStat", "{runnames}"),
        Nanoplot_Dynamic_Histogram_Read_length_html = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_Dynamic_Histogram_Read_length.html"),
        Nanoplot_HistogramReadlength_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_HistogramReadlength.png"),
        Nanoplot_LengthvsQualityScatterPlot_dot_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_LengthvsQualityScatterPlot_dot.png"),
        Nanoplot_LengthvsQualityScatterPlot_kde_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_LengthvsQualityScatterPlot_kde.png"),
        Nanoplot_LogTransformed_HistogramReadlength_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_LogTransformed_HistogramReadlength.png"),
        Nanoplot_NanoPlot_20200515_2130_log = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_NanoPlot_20200515_2130.log"),
        Nanoplotreport_html = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_NanoPlot-report.html"),
        Nanoplot_NanoStats_txt = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_NanoStats.txt"),
        Nanoplot_Weighted_HistogramReadlength_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_Weighted_HistogramReadlength.png"),
        Nanoplot_Weighted_LogTransformed_HistogramReadlength_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_Weighted_LogTransformed_HistogramReadlength.png"),
        Nanoplot_Yield_By_Length_png = os.path.join("QC", "runs", "NanoPlot", "{runnames}", "{runnames}_Yield_By_Length.png")
    run:
        # try:
        #     os.makedirs(os.path.join("QC", "runs", "MinionQC"))
        # except FileExistsError:
        #     pass
        args = {
        "input":input.seq_summary,
        # "outputMin":os.path.join("QC", "runs", "MinionQC"),
        # "minionQCpath":"/media/utlab/DATA_HDD1/Nanopore_metagenomics/Softwares_for_analysis/minion_qc/MinIONQC.R",
        # "minionQCpath":snakemake.config["minionQCpath"]
        "outputNanoS":os.path.join("QC", "runs", "NanoStat"),
        "outputNanoP":os.path.join("QC", "runs", "NanoPlot", wildcards.runnames),
        "name": wildcards.runnames,
        "prefix": wildcards.runnames+"_"
        }
        # shift minionQCpath to config
        #  -s makes small figures suitable for export rather than optimised for screen
        # command = "Rscript {minionQCpath} -i {input} -o {outputMin} -s TRUE"
        # shell(command.format(**args))
        command_nanoS = "NanoStat --summary {input} --outdir {outputNanoS} -n {name} --readtype 1D"
        shell(command_nanoS.format(**args))
        command_nanoP = "NanoPlot --summary {input} --outdir {outputNanoP} -p {prefix} --readtype 1D"
        shell(command_nanoP.format(**args))
#
# include run/s that are/were live basecalled or were only available as fastq? USE qcat
#
# rule demultiplex:
#     input:
#         raw_fastq=rules.basecalling.output.basecalled_dir
#     output:
#         demux_dir=directory(os.path.join(config['ROOT'], "qcat_output/demuxd", "{runnames}")),
#         trimmed_dir=directory(os.path.join(config['ROOT'], "qcat_output/trimmed", "{runnames}"))
#     run:
#         args = {
#         "input":input.raw_fastq,
#         "output.demux_dir":output.demux_dir,
#         "output.trimmed_dir":output.trimmed_dir,
#         "kit":config['barcode_kit']
#         }
#         command = "qcat -fastq {input} --barcode_dir {output.demux_dir} --output {output.trimmed_dir} --trim -k {kit} --detect-middle"
#         command = command.format(**args)
#         shell(print)
#         shell(command)
#
# QC of fastq files
# rule sampleQC:
#     input:
#         "each sample fastq"
#     output:
#         "reports per sample"
#     run:
#         Nanoplot
#         Nanostat
#
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
