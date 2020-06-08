
This is an experimental snakemake pipeline to carry out the entire analysis of my Chronic Pancreatitis (CP) microbiome data and obtained by nanopore sequencing.

> This is not intended to be a ready to use pipeline and I doubt that I will be revisiing it for updates in the future. I learnt Snakemake from scratch through the course of this project and would have found a repo like this immensely useful. That is the main motivation I had to create this repository. Most of the code can be easily understood and reused.<br/>


<p align="center">
  <img src="rulegraph.png" width="500" />
</p>

# Usage
The pipeline includes, data organization, QC, Basecalling (Guppy for nanopore), Demultiplexing, Trimming and filtering, Taxonomic classification (using Kraken2 and Centrifuge).

Be sure to look through the comments in the Snakefile to choose which parts to use and what to edit out. (I may add the info there to README later). Currently it is written to handle a case where one has multiplexed different samples across multiple runs. In such a case, the wildcards change from sample aggregation onwards. So, the workflow will have to be run for the first time with the targets using the sample wildcards commented out. Once that is done, the sample collectSamples rules and onwards can be run. Otherwise, a missing output exception will be thrown. Also, if the upstream run fastq file is edited make sure that the first part of the workflow is redone as needed.


# Tools used

Below is a list of tools used by this pipeline with a link to their page/repo where information about installation and usage can be found.

- [Guppy](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html) (you need to be part of the nanopore community for access to official docs and [download](https://community.nanoporetech.com/downloads))
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [qcat](https://github.com/nanoporetech/qcat)
- [Kraken 2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual)
- [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)
- [NanoFilt](https://github.com/wdecoster/nanofilt/)

If anybody will find this useful if polished, please let me know. Also, if you are a veteran and somehow come across this repo I will be thrilled to hear your sugesstions and more so if someone has the time to help me polish it. Finally, if you find any part of this repo useful for your work, do let me know. I will be very happy to hear!

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
    ├── kraken2_BacArchViFunProt<br/>
    ├   ├── report<br/>
    ├   └── result<br/>
    ├── kraken2_humandb<br/>
    ├   ├── report<br/>
    ├   ├── result<br/>
    ├   └── unclassified<br/>
    └── kraken2_Minidb<br/>
        ├── report<br/>
        └── result<br/>

### NanoFilt added<br/>
input of rule kraken2_human and following need to be edited depending on which fastq files (filtered Q > 7, Q > 10 or unfiltered) are to be used for downstream analysis. Edit this rule (filterSamples) as needed. Refere to the NanoFilt repository for more details.
