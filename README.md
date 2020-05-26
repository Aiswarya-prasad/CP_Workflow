# CP_Workflow
The aim is to create a snakemake workflow to carry out the entire analysis of my Chronic Pancreatitis (CP) microbiome data and compare it with controls.
> This is not intended to be a ready to use pipeline and I doubt that I will be revisiing it for updates in the future. I learnt Snakemake from scratch through the course of this project and would have found a repo like this immensely useful. That is the main motivation I had to create this repository. Most of the code can be easily understood and reused.<\br>

The workflow aims to include, data organization, QC, Basecalling (Guppy for nanopore), Demultiplexing, Trimming and filtering, Taxonomic classification (using Kraken2 and Centrifuge) and comparitive analysis.
It also aims to include a compilation of all the analysis in the form of my Master thesis written in Latex.

directory structure of output as a classification<\br>
classified<\br>
└── <sample#><\br>
    ├── bracken<\br>
    │   ├── genus_report<\br>
    │   └── species_report<\br>
    ├── centrifuge<\br>
    │   ├── report<\br>
    │   └── result<\br>
    ├── kraken2_BacArchViFunProt<\br>
    │   ├── report<\br>
    │   └── result<\br>
    ├── kraken2_humandb<\br>
    │   ├── report<\br>
    │   ├── result<\br>
    │   └── unclassified<\br>
    └── kraken2_Minidb<\br>
        ├── report<\br>
        └── result<\br>
