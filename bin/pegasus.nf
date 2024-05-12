#!/usr/bin/env nextflow
// initialzie and declare input channels
nextflow.enable.dsl=2
params.publish_dir = './Results'

params.nanoPath = ""
params.shortPath1 = ""
params.shortPath2 = ""
params.buscoPath = ""
params.quastPathFNA = ""
params.quastPathGFF = ""
params.phvDatabase = ""
params.ref1 = ""

params.centrifugeRscript = "/users/j/l/jlubkowi/scratch/PEGASUS/run_centrifuge_clean.R"
centrifugeRscript = Channel.fromPath( params.centrifugeRscript, checkIfExists: true)

phvDatabase = Channel.fromPath( params.phvDatabase, checkIfExists: true)
shortreads1_ch = Channel.fromPath( params.shortPath1, checkIfExists: true )
shortreads2_ch = Channel.fromPath( params.shortPath2, checkIfExists: true )
longread_ch = Channel.fromPath( params.nanoPath , checkIfExists: true)
quastpathGFF_ch = Channel.fromPath( params.quastPathGFF, checkIfExists: true)
quastpathFNA_ch = Channel.fromPath( params.quastPathFNA, checkIfExists: true)
buscopath_ch = Channel.fromPath( params.buscoPath, checkIfExists: true)
ref1_ch = Channel.fromPath( params.ref1, checkIfExists: true)

// pipeline processes

process NANOPLOT {
    publishDir "${params.publish_dir}", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.0--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.0--pyhdfd78af_0' }"   

    input:
        path nano_read 
    output:
        path 'Nanoplot' 
    script:
        """
        mkdir Nanoplot
        cp ${nano_read} Nanoplot
        cd Nanoplot
        NanoPlot -t 40 --fastq ${nano_read}
        """
}

process FASTQC {
    publishDir "${params.publish_dir}", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'biocontainers/fastqc:0.11.9--0' }"

    input:
        path short_reads1 
        path short_reads2
        path nano_reads 

    output:
        path 'FastQC'
    script: 
    """
        mkdir FastQC
        fastqc --noextract --nogroup -o FastQC ${short_reads1}/*
        fastqc --noextract --nogroup -o FastQC ${short_reads2}/*
        fastqc --noextract -o FastQC ${nano_reads} -t 8
    """
}

process FASTP_LONG {
     publishDir "${params.publish_dir}/Fastp_long", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_3':
        'biocontainers/fastp:0.23.4--hadf994f_3' }"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        path long_reads

    output:
        path 'long_reads_clean.fastq.gz'

    script:
    """
    fastp -i ${long_reads} -o long_reads_clean.fastq.gz --thread 10
    """
}

process LONG_CENTRIFUGE {
    publishDir "${params.publish_dir}/long_centrifuge", mode: 'copy'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
        path longreads 
        path database 

    output:
        path 'long_cent_out.txt'
        path 'longreads.fastq'

    script:
    """
        cp ${longreads} \$PWD/longreads.fastq.gz
        cp ${database}/* \$PWD

        centrifuge -q -x p+h+v -U longreads.fastq.gz --threads 8 -S long_cent_out.txt --report-file long_cent_out.tsv --min-hitlen 250
            
        gunzip longreads.fastq.gz
    """
}

process R_PROCESSING_LONG {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'biocontainers/r-tidyverse%3A1.2.1' }"

    input:
        path rscript
        path long_input_txt
        path long_fasta
    output:
        path 'not_contam.txt'
        path 'longreads.fastq'
    script:
    """
    cp ${rscript} \$PWD/centrifugeClean.R

    Rscript --vanilla centrifugeClean.R long_cent_out.txt longreads
    """
}

process LONG_SEQTK {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"
        
    input:
        path sequence
        path contam 
    output:
        path 'clean_longreads.fastq.gz'
    script:
    """
    seqtk subseq longreads.fastq not_contam.txt > clean_longreads.fastq

    gzip clean_longreads.fastq
    """
}

process FLYE {

    publishDir "${params.publish_dir}/Flye", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
        path nano_read 

    output:
	    file 'flye/assembly.fasta'
        file 'contig_stats.txt'
    script:
        """
        flye --nano-hq ${nano_read} --out-dir flye --threads 40
        cp \$PWD/flye/assembly.fasta \$PWD/assembly.fasta
        tail \$PWD/flye/flye.log > contig_stats.txt
        """
}

process FASTP_SHORT {
     publishDir "${params.publish_dir}", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_3':
        'biocontainers/fastp:0.23.4--hadf994f_3' }"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        path shortreads1 
        path shortreads2 
    output:
        path 'trimmed_short_reads'

    script:
    """
    mkdir trimmed_short_reads
    cp -r ${shortreads1} trimmed_short_reads
    cp -r ${shortreads2} trimmed_short_reads
    cd trimmed_short_reads
    fastp -i ${shortreads1}/*_R1_001.fastq.gz -I ${shortreads1}/*_R2_001.fastq.gz -o short_reads.1_1.trim.fastq.gz -O short_reads.1_2.trim.fastq.gz --thread 40
    fastp -i ${shortreads2}/*_R1_001.fastq.gz -I ${shortreads2}/*_R2_001.fastq.gz -o short_reads.2_1.trim.fastq.gz -O short_reads.2_2.trim.fastq.gz --thread 40
    """
}

process SHORT_CENTRIFUGE {
    publishDir "${params.publish_dir}/short_centrifuge", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
        path short_reads 
        path database 

    output:
        path 'short_cent'
        path 'short_cent_out.txt'
    
    script:
        """
        cat ${short_reads}/*.1_1.trim.fastq.gz ${short_reads}/*.2_1.trim.fastq.gz > short_reads_1.trim.fastq.gz
        cat	${short_reads}/*.1_2.trim.fastq.gz ${short_reads}/*.2_2.trim.fastq.gz > short_reads_2.trim.fastq.gz
        cp ${database}/* \$PWD
        centrifuge -q -x p+h+v -1 short_reads_1.trim.fastq.gz -2 short_reads_2.trim.fastq.gz --threads 8 -S short_cent_out.txt --report-file short_cent_out.tsv --min-hitlen 20

        gunzip short_reads_1.trim.fastq.gz

        gunzip short_reads_2.trim.fastq.gz

        mkdir short_cent
        mv short_reads_1.trim.fastq short_cent
        mv short_reads_2.trim.fastq short_cent
        """
}

process R_PROCESSING_SHORT {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'biocontainers/r-tidyverse%3A1.2.1' }"

    input:
        path rscript
        path short_inputs 
        path short_cent 
    output:
        path 'not_contam.txt'
        path 'short_cent'
    script:
    """
    cp ${rscript} \$PWD/centrifugeClean.R

    Rscript --vanilla centrifugeClean.R short_cent_out.txt short_reads
    """
}

process SHORT_SEQTK {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
        path sequences 
        path not_contam 

    output:
        path 'short_centrifuge_results'
    script:
    """
    seqtk subseq short_cent/short_reads_1.trim.fastq not_contam.txt > clean_short_reads_1.trim.fastq
    gzip clean_short_reads_1.trim.fastq

    seqtk subseq short_cent/short_reads_2.trim.fastq not_contam.txt > clean_short_reads_2.trim.fastq
    gzip clean_short_reads_2.trim.fastq


    mkdir short_centrifuge_results
    mv clean_short_reads_2.trim.fastq.gz short_centrifuge_results
    mv clean_short_reads_1.trim.fastq.gz short_centrifuge_results
    """
}

process HAPOG {
    publishDir "${params.publish_dir}/Hapog1", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hapog:1.3.7--py39hf95cd2a_0' :
        'biocontainers/hapog:1.3.7--py39hf95cd2a_0' }"

    input:
        path flye 
        path flye_extra 
        path centrifuge_out 

    output:
        path 'hapog.fasta' 

    script: 
    """
    hapog --genome assembly.fasta --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    cp \$PWD/polishing/hapog_results/hapog.fasta \$PWD/hapog.fasta
    """
}


process CHOPPER {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.7.0--hdcf5f25_0' :
        'biocontainers/chopper:0.7.0--hdcf5f25_0' }"

    input:
        path long_reads
    output:
        path 'filtered_reads.fastq.gz'
    script: 
    """
    gunzip -c ${long_reads} | chopper -q 20 -l 500 | gzip > filtered_reads.fastq.gz
    """
}

process MINIMAP {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.26--he4a0461_2' :
        'biocontainers/minimap2:2.26--he4a0461_2' }"
    input:

        path genome 
        path long_reads 
    
    output:
        path 'aligned_long_reads.sam'

    script: 
    """
    minimap2 -ax map-ont --secondary=no ${genome} ${long_reads} > aligned_long_reads.sam
    """
}

process RACON1 {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.5.0--h7ff8a90_1' :
        'biocontainers/racon:1.5.0--h7ff8a90_1' }"
    input:
        path genome 
        path long_reads
        path long_reads_sam

    
    output:
        path 'genome_SLpolished.fasta'

    script: 
    """
    racon ${long_reads} ${long_reads_sam} ${genome} -t 10 > genome_SLpolished.fasta
    """
}


process NTLINK {
    tag "ntLinking"
    publishDir "${params.publish_dir}/NTLink", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ntlink:1.3.9--py39hd65a603_2' :
        'biocontainers/ntlink:1.3.9--py39hd65a603_2' }"

    input:
        path nano_reads 
        path genome 
    output:
        path 'hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa'

    script: 

    """
    cp ${genome} \$PWD/hapog_result.fasta
    cp ${nano_reads} \$PWD/nanopore_raw.fastq.gz
    ntLink_rounds run_rounds target=hapog_result.fasta reads=nanopore_raw.fastq.gz k=32 w=250 t=16 overlap=True rounds=5
    """
}

process HAPOG2 {
    tag "HAPOGing the NTLINK"
    publishDir "${params.publish_dir}/Hapog2", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hapog:1.3.7--py39hf95cd2a_0' :
        'biocontainers/hapog:1.3.7--py39hf95cd2a_0' }"

    input:
        path ntlink 
        path centrifuge_out

    output:
        path 'hapog.fasta' 

    script: 
    """
    hapog --genome hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    mv \$PWD/polishing/hapog_results/hapog.fasta \$PWD/hapog.fasta
    """
}

process SEQTK_RENAME {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
        path genome
    output:
        path 'genome_renamed.fasta'
    script:
    """
    seqtk rename ${genome} n > genome_renamed.fasta
    """
}

process MINIMAP2 {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.26--he4a0461_2' :
        'biocontainers/minimap2:2.26--he4a0461_2' }"
    input:

        path genome 
        path long_reads 
    
    output:
        path 'aligned_long_reads.sam' 
    script: 
    """
    awk -F ' ' '/^>/{print \$1; next} 1' ${genome} > genome_fixed_headers.fasta
    minimap2 -ax map-ont -N0 genome_fixed_headers.fasta ${long_reads} > aligned_long_reads.sam
    """
}

process RACON2 {
    publishDir "${params.publish_dir}/Racon2", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.5.0--h7ff8a90_1' :
        'biocontainers/racon:1.5.0--h7ff8a90_1' }"
    input:
        path genome 
        path long_reads
        path long_reads_sam

    output:
        path 'genome_SL2polished.fasta'

    script: 
    """
    racon ${long_reads} ${long_reads_sam} ${genome} -t 10 > genome_SL2polished.fasta
    """
}

process RAGTAG {
    tag "Ragging and Tagging "
    publishDir "${params.publish_dir}", mode: 'copy' 
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0' :
        'biocontainers/ragtag:2.1.0--pyhb7b1952_0' }"

    input:
        path genome
        path ref1
    output:
        path 'Ragtag'
    script: 
    """
    cp ${ref1} \$PWD/ref1.fasta
    cp ${genome} \$PWD/results.fasta

    ragtag.py scaffold -o Ragtag ref1.fasta results.fasta
    """
}

process BUSCOS {
    tag "Busco"
    publishDir "${params.publish_dir}/Busco", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.5.0--pyhdfd78af_0' :
        'biocontainers/busco:5.5.0--pyhdfd78af_0' }"
    input:
        path final_genome_path 
        path ragtag_genome
        path busco 
    output:
        path 'busco_genome'
        path 'busco_ragtag'
    script:
    """
    cp ${final_genome_path} \$PWD/results.fasta
    cp ${ragtag_genome}/ragtag.scaffold.fasta \$PWD/ragtag.fasta
    
    busco -m genome -i results.fasta -o busco_genome -l ${busco}

    busco -m genome -i ragtag.fasta -o busco_ragtag -l ${busco}
    """
}

process QUASTGENOME {
    tag "Quasting"

    publishDir "${params.publish_dir}", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path final_genome_path 
        path quastGenomeGFF 
        path quastGenomeFNA 
    output: 
        path 'quast_results'
    script:
    """
    cp ${quastGenomeGFF} \$PWD/input.gff
    cp ${quastGenomeFNA} \$PWD/input.fna
    cp ${final_genome_path} \$PWD/genome.fasta
    
    quast genome.fasta -r input.fna -g input.gff --large -t 10
    """
}

process QUASTRAGTAG {
    tag "Quasting"

    publishDir "${params.publish_dir}", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path ragtag_genome 
        path quastGenomeGFF 
        path quastGenomeFNA

    output: 
        path 'ragtag_quast_results'
    script:
    """
    cp ${quastGenomeGFF} \$PWD/input.gff
    cp ${quastGenomeFNA} \$PWD/input.fna
    cp ${ragtag_genome}/ragtag.scaffold.fasta \$PWD/ragtag.fasta
        
    quast ragtag.fasta -r input.fna -g input.gff --large -t 10

    cp -r quast_results ragtag_quast_results
    """
}
 
process BAMCONVERSION {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
        path genome 
        path centrifuge_out 
    
    output:
        path 'short_reads_mapped2_sorted.bam' 
        path 'short_reads_mapped2_sorted.bam.bai'

    script: 
    """
    # index
    bwa-mem2 index -p mem2 ${genome}

    # map
    bwa-mem2 mem mem2 -t 40 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    #convert to bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    #sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam

    #index
    samtools index short_reads_mapped2_sorted.bam
    """
}

process MOSDEPTH {
    tag "Mosdepthing"
    publishDir "${params.publish_dir}/Mosdepth", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.6--hd299d5a_0' :
        'biocontainers/mosdepth:0.3.6--hd299d5a_0' }"

    input:
        path genome 
        path index 
    output:
        path 'mosdepth_result.*'
    script: 
    
    """
    mosdepth -t 8 mosdepth_result ${genome}
    """
}

process MULTIQC {
    publishDir "${params.publish_dir}/MultiQC", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"
    
    input:
        path pore_reads
        path pore_data
        path nano_ch 
        path fastqc_ch
        path busco_rag
        path busco_gen
        path quast_gen_ch 
        path quast_rag_ch
        path mosdepth_txt
 

    output:
        path 'multiqc_report.html'
    script:
        """
        multiqc . -d -dd 3
        """
}

workflow {
    // Stats
    nano_ch = NANOPLOT(longread_ch)
    fastqc_ch = FASTQC(shortreads1_ch, shortreads2_ch, longread_ch)

    // Short Read prep
    fastp_short_ch = FASTP_SHORT(shortreads1_ch, shortreads2_ch)
    SHORT_CENTRIFUGE(FASTP_SHORT.out, phvDatabase)
    R_PROCESSING_SHORT(centrifugeRscript, SHORT_CENTRIFUGE.out)
    SHORT_SEQTK(R_PROCESSING_SHORT.out)

    // Long Read prep
    fastp_long_ch = FASTP_LONG(longread_ch)
    LONG_CENTRIFUGE(FASTP_LONG.out, phvDatabase)
    R_PROCESSING_LONG(centrifugeRscript, LONG_CENTRIFUGE.out)
    LONG_SEQTK(R_PROCESSING_LONG.out)
    CHOPPER(LONG_SEQTK.out)

    // Assembly
    FLYE(LONG_SEQTK.out)

    HAPOG(FLYE.out, SHORT_SEQTK.out)    
    MINIMAP(HAPOG.out, CHOPPER.out)
    RACON1(HAPOG.out, CHOPPER.out, MINIMAP.out)

    NTLINK(longread_ch, RACON1.out)

    HAPOG2(NTLINK.out, SHORT_SEQTK.out)
    SEQTK_RENAME(HAPOG2.out)
    MINIMAP2(SEQTK_RENAME.out, CHOPPER.out)
    RACON2(SEQTK_RENAME.out, CHOPPER.out, MINIMAP2.out)

    RAGTAG(RACON2.out, ref1_ch)

    // Quality Checks
    busco_ch = BUSCOS(RACON2.out, RAGTAG.out, buscopath_ch)
    quast_gen_ch = QUASTGENOME(RACON2.out, quastpathGFF_ch, quastpathFNA_ch)
    quast_rag_ch = QUASTRAGTAG(RAGTAG.out, quastpathGFF_ch, quastpathFNA_ch)
    BAMCONVERSION(RACON2.out, SHORT_SEQTK.out)
    mosdepth_ch = MOSDEPTH(BAMCONVERSION.out)
    MULTIQC(fastp_long_ch, fastp_short_ch, nano_ch, fastqc_ch, busco_ch, quast_gen_ch, quast_rag_ch,mosdepth_ch)
}