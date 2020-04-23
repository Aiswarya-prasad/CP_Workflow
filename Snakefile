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


# --- Some rules --- #

rule all:
    input:
        expand(os.path.join("fastq", "{runnames}.fastq"), runnames=config['runnames'])
    threads: 2

# JUST 2 RULES THAT WORK
# rule all:
#     input:
#         expand(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"), runnames=RUNNAMES)
#     threads:2
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
#
#
rule basecalling:
    input:
        raw_dir=os.path.join(config['RAWDIR'], "{runnames}")
    output:
        # basecalled_dir=directory(os.path.join(config['RAWDIR'], "guppy_output", "{runnames}"))
        run_fastq=os.path.join("fastq", "{runnames}.fastq")
    run:
        # if recursive is enabled, I can give the run dir which has sub runs..(Expt 4) (--min_qscore 7 is already default)
        # conditionally allow for --resume figure out how later
        # there is a recent issue April2020 with barcode trimming in guppy_basecaller so use qcat for demult
        # dna_r9.4.1_450bps_hac.cgf for FLO-MIN106 and SQK-LSK109 combination
        guppy_output_dir = os.path.join(config['RAWDIR'], "guppy_output", wildcards.runnames)
        try:
            os.makedirs(guppy_output_dir)
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
           try:
               shell("rename "+guppy_output_dir+"/pass/*.fastq "+guppy_output_dir+"/{runnames}.fastq")
               shell("rsync -v "+guppy_output_dir+"/pass/"+wildcards.runnames+".fastq fastq/"+wildcards.runnames+".fastq")
           except:
               print("no basecalling happened")
               shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
           else:
               pass
           # shell("rename "+guppy_output_dir"+/pass/*.fastq "+guppy_output_dir"+/{runnames}.fastq")
           # shell("rsync -v "+guppy_output_dir"+/pass/*.fastq fastq")
        else:
            print("No log file to resume from. Starting fresh instance of basecallig")
            command = "guppy_basecaller --input_path {input} --save_path {output_dir} --flowcell FLO-MIN106 --kit SQK-LSK109 --recursive --records_per_fastq 0 --calib_detect --qscore_filtering"
            command = command.format(**args)
            try:
                shell(command)
                shell("rename "+guppy_output_dir+"/pass/*.fastq "+guppy_output_dir+"/{runnames}.fastq")
                shell("rsync -v "+guppy_output_dir+"/pass/"+wildcards.runnames+".fastq fastq/"+wildcards.runnames+".fastq")
            except:
                print("no basecalling happened")
                shell("touch "+"fastq"+"/"+wildcards.runnames+".fastq")
            else:
                pass

#
#
# For run1 copy seq summary into output dir
# find sequencing summary
# rule runQC:
#     input:
#         seq_summary=os.path.join(rules.basecalling.output.basecalled_dir, "sequencing_summary.txt")
#     output:
#         "reports per run"
#     shell:
#         MinionQC
#         Nanocomp
#
# include run/s that are/were live basecalled or were only available as fastq? USE qcat
#
# rule demultiplex:
#     input:
#         raw_fastq=rules.basecalling.output.basecalled_dir
#     output:
#         demux_dir=directory(os.path.join(config['RAWDIR'], "qcat_output/demuxd", "{runnames}")),
#         trimmed_dir=directory(os.path.join(config['RAWDIR'], "qcat_output/trimmed", "{runnames}"))
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
