![VBRN Logo](./logo_landscape.png)

# Hybrid Genome Assembly 
# HyGBuild

### Authors
[Emily Curd](https://scholar.google.com/citations?user=uGHWHbgAAAAJ&hl=en&oi=ao) and James (Jax) Lubkowitz

### Affiliations
[Vermont Biomedical Research Network](https://vbrn.org), [Vermont Integrated Genomics Resource](https://www.med.uvm.edu/vigr/home), [University of Vermont](https://www.uvm.edu), [St. Lawrence University](https://www.stlawu.edu)

### Funding
Research reported in this repository was supported by an Institutional Development Award (IDeA) from the National Institute of General Medical Sciences of the National Institutes of Health under grant number P20GM103449. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of NIGMS or NIH.

### Installation and Depedencies
This pipeline requires Singularity, Nextflow and Hapo-G.
- [Nextflow Installation](https://www.nextflow.io/docs/latest/getstarted.html)

- Hapo-g can be pulled with the rest of the pipeline or manually via [cloning the Github](https://github.com/institut-de-genomique/HAPO-G)
    - Note: Hapo-G relies on bwa and samtools 

```
git clone https://github.com/institut-de-genomique/HAPO-G hapog
cd hapog
bash build.sh

# check bwa and samtools installation 
    which bwa
    which samtools
```

Once completed, clone the pipeline repository:
```
git clone https://github.com/jaxlub/TODDO
```

look into adding it to PATH/Executable

### Using ____




### Pipeline Contents
![Pipeline Diagram](./pipeline_diagram.png)
Color corresponds to process that can be run in parralel




### Parameters

- nanoPath: path to directory containing nanopore longreads
- shortPath1
- shortPath2
- trimpath1: path to directory containing first batch of shortreads
- trimpath2: path to directory containing second batch of shortreads
- hapog: path to hapog executable
- quastPathGFF: path to GFF genome for comparison 
- quastPathFNA: path to FNA genome for comparison 
- buscoPath: Path to reference genome for Busco comparison
- ref1: Path to reference genome for ragtag assembly

### Example usage
```
nextflow run  /users/j/l/jlubkowi/scratch/bullhead_project/HL4/main.nf -resume \
            --nanoPath '/users/j/l/jlubkowi/scratch/bullhead_project/HL4/reads/long_reads/*' \
            --shortPath1 '/users/j/l/jlubkowi/scratch/bullhead_project/HL4/reads/short_reads1/*' \
            --shortPath2 '/users/j/l/jlubkowi/scratch/bullhead_project/HL4/reads/short_reads2/*' \
            --trimpath1 '/users/j/l/jlubkowi/scratch/bullhead_project/HL4/reads/short_reads1' \
            --trimpath2 '/users/j/l/jlubkowi/scratch/bullhead_project/HL4/reads/short_reads2/' \
            --hapog '/users/j/l/jlubkowi/scratch/bullhead_project/hapog/hapog.py' \
            --quastPathGFF '/users/j/l/jlubkowi/scratch/bullhead_project/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.gff' \
            --quastPathFNA '/users/j/l/jlubkowi/scratch/bullhead_project/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.fna' \
            --buscoPath '/gpfs1/home/e/g/eguswa/scratch/bullhead/actinopterygii_odb10' \
            --ref1 '/users/j/l/jlubkowi/scratch/bullhead_project/Ameiurus_melas_genome/GCA_012411365.1_AMELA_1.0_genomic.fna' 
```

### Output 
Pipeline outputs all notable results to a Results directory. 
All intermediate information and files are stored in the nextflow work directories and now published. 

Published Results Directory Structure
```
Results
    -- Busco
        -- busco_genome
            -- Contains Busco results of hapog-built genome
        -- busco_ragtag
            -- Contains Busco results of ragtag-scaffolded genome
    -- FastQC
        -- Contains FastQC zip and html files
    -- Flye
        -- flye
            -- Contains assembly.fasta, an intermediate genome
    -- Hapog1
        -- polishing
            -- hapog
                -- Contains hapog.fasta an intermediate genome
    -- Hapog2
        -- polishing
            -- hapog
                -- Contains hapog.fasta an intermediate genome
    -- long_centrifuge
        -- Contains cleaned longreads in fastq.gz
    -- Mosdepth
        -- Contains mosdepth result files
    -- Nanoplot
        -- Contains Nanoplot png and html files
    -- NTLink
        -- Contains hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa, an intermediate genome
    -- RagTag
        -- ragtag_results
            -- Contains all ragtag results
    -- short_centrifuge
        -- short_centrifuge_results
            -- Contains cleaned shortreads in fastq.gz
    -- Trimmomatic
        -- trimmomatic
            -- Contains trimmed and untrimmed shortreads
``` 


### Refrences and Tool Citations

#### Tools
##### Busco
Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323

##### Centrifuge
Kim D, Song L, Breitwieser FP, and Salzberg SL. Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome Research 2016

http://ccb.jhu.edu/software/centrifuge/

##### FastQC
Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Data accessed: 12/6/2023



Flye
##### Hapo-G
Jean-Marc Aury, Benjamin Istace, Hapo-G, haplotype-aware polishing of genome assemblies with accurate reads, NAR Genomics and Bioinformatics, Volume 3, Issue 2, June 2021, lqab034, https://doi.org/10.1093/nargab/lqab034

MultiQC
##### Ntlink
Coombe L, Li JX, Lo T, Wong J, Nikolic V, Warren RL and Birol I. LongStitch: High-quality genome assembly correction and scaffolding using long reads. bioRxiv. 2021;2021.06.17.448848. doi: https://doi.org/10.1101/2021.06.17.448848

##### Trimmomatic
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
Published online 2014 Apr 1. doi: 10.1093/bioinformatics/btu170

##### Ragtag

##### Quast

#### Containers
- Mosdepth v0.24.0 Container \
Handler, Dominik. 'dominik-handler/AP_singu:mosdepth'. 
Data Accessed: 11/29/2023
Commit: e55ba3e87926f1ae97d6e85fa8284f37f594a833
https://singularityhub.github.io/singularityhub-archive/containers/dominik-handler-AP_singu-mosdepth/

