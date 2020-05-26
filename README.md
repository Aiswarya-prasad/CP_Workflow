# CP_Workflow
The aim is to create a snakemake workflow to carry out the entire analysis of my Chronic Pancreatitis (CP) microbiome data and obtained by nanopore sequencing.

> This is not intended to be a ready to use pipeline and I doubt that I will be revisiing it for updates in the future. I learnt Snakemake from scratch through the course of this project and would have found a repo like this immensely useful. That is the main motivation I had to create this repository. Most of the code can be easily understood and reused.<br/>

The pipeline includes, data organization, QC, Basecalling (Guppy for nanopore), Demultiplexing, Trimming and filtering, Taxonomic classification (using Kraken2 and Centrifuge) and a seperate script for data analysis which can be called by a different rule.

If anybody will find this useful if polished, please let me know. Also, if you are a veteran and somehow come across this repo I will be thrilled to hear your sugesstions and more so if someone has the time to help me polish it. Finally, if you find any part of this repo useful for your work, do let me know. I will be very happy to hear!

directory structure of output as a classification<br/>
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
