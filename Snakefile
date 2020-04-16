# Main Workflow - CP project analysis
#
# Contributors: @aiswaryaprasad

# --- Importing Configuration Files --- #
configfile: "config.yaml"



# --- Some rules --- #

rule parse_metadate:
    # somehow parse metadate and make it available to all the other rules as needed

rule find_fast5:
    input:
        "path to Run dir"
    output:
        "path to fast5 files"
    script:
        "scripts/find_fast5.py"

rule raw_QC:
    input:
        "path to sequencing_summary.txt"
    output:
        "figures",
        "text"
    shell:
        "program for doing that minionQC (?????)"

        ######## --- will depend on --- ########
        rule find_sumary:
            "path to Run dir"
        output:
            "path to seq summary"
        script:
            "scripts/find_sumary.py"

rule basecalling:
    input:
        "the fast5 files"
    output:
        "consider syntax and documentation"
    shell:
        "write the guppy_basecaller"

rule fastq_QC_run:
    input:
        "basecalled fastq files per run"
    output:
        "figures",
        "text or csv"
    shell:
        "appropriate program to do QC. If minionQC works, drop this"

rule demultiplex:
    input:
        "basecalled fastq files per run"
    output:
        "all barcodes for each run"
    shell:
        "guppy_barcoder (ok? -- other options)"

        ######## --- will depend on --- ########
        # sample, run, barcode correlation from config or metadata (latter better)
        rule seggregate_fastq:
            input:
                "result of demux per run"
            output:
                "fastq for each sample"
            script:
                "script/seggregate_fastq.py"

rule run_kraken2:
    input:
        "fatq files per sample"
    output:
        "reports",
        "results"
    shell:
        "call run_kraken2"

rule run_centrifuge:
    input:
        "fatq files per sample"
    output:
        "reports",
        "results"
    shell:
        "call run_centrifuge"

###########--figure out how to do comparative analysis--###########
