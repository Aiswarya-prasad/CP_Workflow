
This is an experimental snakemake pipeline to carry out the entire analysis of my Chronic Pancreatitis (CP) microbiome data and obtained by nanopore sequencing.

> This is not (yet) intended to be a ready-to-use pipeline. I learnt Snakemake from scratch through the course of this project and would have found a repo like this immensely useful. That is the main motivation I had to create this repository. Most of the code can be easily understood and reused.<br/>


<p align="center">
  <img src="rulegraph.png" width="500" />
</p>


# Tools used

Below is a list of tools used by this pipeline with a link to their page/repo where information about installation and usage can be found.

- [Guppy](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html) (you need to be part of the nanopore community for access to official docs and [download](https://community.nanoporetech.com/downloads))
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [qcat](https://github.com/nanoporetech/qcat)
- [NanoFilt](https://github.com/wdecoster/nanofilt/)
- [Kraken 2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual)
- [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)

If anybody will find this useful if polished, please let me know. Also, if you are a veteran and somehow come across this repo I will be thrilled to hear your sugesstions and more so if someone has the time to help me polish it. Finally, if you find any part of this repo useful for your work, do let me know. I will be very happy to hear!

# Setup and use</br>

Firstly, all these tools need to be installed and working. Some of this are possible to install through conda. Snakemake can handle this and I will soon add the conda environment set up for this. For now, they need to be installed and set up manually.

In the config.yaml file, make sure to configure the following:

-RAWDIR: Path to directory containing links to fast5 directory (or) path to fast5 directory to be used by Guppy
-ROOT: Path to project directory to be used by shell commands and scripts if necessary
runnames: List of run names eg. ['Run0', 'Run1', 'Run2', 'Run3', 'Run4']
-barcode_kit: 'NBD103/NBD104' parameter used by qcat
samples: List of sample IDs eg. ['01', '03', '06', '07', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22']
-sample_dict: a nested dictionary mapping runs to barcodes and sample IDs. ie. the first level of keys are run names, within each dict which is a value for the key run name, is a dictonary where keys are barcode numbers and values are sample ID associated with each barcode in that particular run. *This is useful for demultiplexing each run and separating out the files based on sample ID rather than barcode as barcodes may not be unique unless coupled with run names leading to long messy names.* **leave this out and edit Snakefile accordingly if this is not desired**
eg.
sample_dict:
    'Run0':
      '04': '01'
    'Run1':
      '02': '03'
      '03': '13'
    'Run2':
      '04': '06'
      '05': '07'
      '06': '10'
      '07': '14'
    'Run3':
      '08': '11'
      '09': '12'
      '10': '15'
      '11': '17'
      '12': '18'
    'Run4':
      '01': '19'
      '02': '20'
      '03': '21'
      '04': '22'
-kraken_db: '$MINIKRAKEN_DB' environment varaible if it is configured in ~/.bashrc. If not, write the path to the kraken database here.
-centrifuge_db: path to the centrifuge database


make sure that these programs are configured properly and added to the path variable in ~/.bashrc as is appropriate.

--------------------------------------------------------------------------------------------------------------------------

Most of the information about output and input formats can be found as comments in the code.

### directory structure of output after classification<br/>
classified<br/>
└── <sample#><br/>
    ├── bracken<br/>
    ├   ├── genus_report<br/>
    ├   └── species_report<br/>
    ├── centrifuge<br/>
    ├   ├── report<br/>
    ├   └── result<br/>
    └── kraken2_Minidb<br/>
        ├── report<br/>
        └── result<br/>

### NanoFilt added<br/>
input to rule running kraken2 and centrifuge can be edited depending on which fastq files (filtered Q > 7, Q > 10 or unfiltered) are to be used for downstream analysis. Edit this rule (filterSamples) as needed. Refere to the NanoFilt repository for more details.
