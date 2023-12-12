#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    tag "FASTQCing"
    publishDir "${params.publish_dir}/FastQC", mode: 'copy'

    conda "bioconda::fastqc=0.11.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'biocontainers/fastqc:0.11.9--0' }"

    input:
        path short_reads1 // channel of first batch of short read files PATH/TO/SHORT/READS/*
        path short_reads2 // channel of first batch of short read files PATH/TO/SHORT/READS/*
        path nano_reads // channel of long reads files PATH/TO/LONG/READS/*

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

process LONG_CENTRIFUGE {
    tag "Centrifuging Long reads"
    publishDir "${params.publish_dir}/long_centrifuge", mode: 'copy'
    
    input:
        path rscript // path to the r-Script
        path longreads // path for longread
        path database // Path to p+h+v indexes
        path rimage // Path to R Image 
        path centrifugeimage // Path to centrifuge image
        path seqtkimage // Path to seqtkimage

    output:
        path 'clean_longreads.fastq.gz'

    script:
    """
        pwds=\$(echo \$PWD)
        cp ${longreads} \$pwds/longreads.fastq.gz
        cp ${database}/* \$pwds
        cp ${rscript} \$pwds/centrigugeClean.R
        
        module load singularity
        singularity exec --bind \$pwds ${centrifugeimage} centrifuge -q -x p+h+v -U longreads.fastq.gz --threads 8 -S ${longreads}cent_out.txt --report-file ${longreads}cent_out.tsv --min-hitlen 250
        
        singularity exec --bind \$pwds ${rimage} Rscript --vanilla centrigugeClean.R ${longreads}cent_out.txt longreads
        
        gunzip longreads.fastq.gz

        file_name=\$(echo '${longreads}' | cut -f 1-2 -d '.')


        singularity exec --bind \$pwds ${seqtkimage} seqtk subseq longreads.fastq not_contam.txt > clean_longreads.fastq

        gzip clean_longreads.fastq
    """
}

process SHORT_CENTRIFUGE {
    tag "Centrifuging"
    publishDir "${params.publish_dir}/short_centrifuge", mode: 'copy'

    input:
        path rscript // path to Centrifuge r-Script
        path trimmomatic // path to trimmomatic result
        path database // Path to p+h+v.tar.gz
        path rimage // Path to R Image 
        path centrifugeimage // Path to centrifuge image
        path seqtkimage // Path to seqtkimage

    output:
        path 'short_centrifuge_results'
    
    script:
        """
        module load singularity
        cat ${trimmomatic}/*.1_1.trim.fastq.gz trimmomatic/*.2_1.trim.fastq.gz > short_reads_1.trim.fastq.gz
        cat	${trimmomatic}/*.1_2.trim.fastq.gz trimmomatic/*.2_2.trim.fastq.gz > short_reads_2.trim.fastq.gz
        pwds=\$(echo \$PWD)
        cp ./PHVindexes/* \$pwds
        cp ${rscript} \$pwds/centrifugeClean.R

        singularity exec --bind \$pwds ${centrifugeimage} centrifuge -q -x p+h+v -1 short_reads_1.trim.fastq.gz -2 short_reads_2.trim.fastq.gz --threads 8 -S short_reads_out.txt --report-file short_reads_out.tsv --min-hitlen 20

        singularity exec --bind \$pwds ${rimage} Rscript --vanilla centrifugeClean.R short_reads_out.txt short_reads

        gunzip short_reads_1.trim.fastq.gz

        singularity exec --bind \$pwds ${seqtkimage} seqtk subseq short_reads_1.trim.fastq not_contam.txt > clean_short_reads_1.trim.fastq

        gzip clean_short_reads_1.trim.fastq

        gunzip short_reads_2.trim.fastq.gz

        singularity exec --bind \$pwds ${seqtkimage} seqtk subseq short_reads_2.trim.fastq not_contam.txt > clean_short_reads_2.trim.fastq

        gzip clean_short_reads_2.trim.fastq

        mkdir short_centrifuge_results
        mv clean_short_reads_2.trim.fastq.gz short_centrifuge_results
        mv clean_short_reads_1.trim.fastq.gz short_centrifuge_results
        """
}


process NANOPLOT {
    tag "Nanoplotting"

    publishDir "${params.publish_dir}/Nanoplot", mode: 'copy'


    conda "bioconda::nanoplot=1.41.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.0--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.0--pyhdfd78af_0' }"   

    input:
        path nano_read // channel of long reads files PATH/TO/LONG/READS/*
    output:
        path 'nanoplot' 
    script:
        """
        mkdir nanoplot
        cp ${nano_read} nanoplot
        cd nanoplot
        NanoPlot -t 40 --fastq ${nano_read}
        """
}


process FLYE {
    tag "Flyeing"

    // Only published result from proccess is assembly.fasta found in workingDir/flye/assembly.fasta and moved to ./results/Flye
    publishDir "${params.publish_dir}/Flye", mode: 'copy'

    conda "bioconda::flye=2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
        path nano_read // channel of long reads files PATH/TO/LONG/READS/*

    output:
	    path 'flye/assembly.fasta'
    script:
        """
        flye --nano-hq ${nano_read} --out-dir flye --threads 40
        """
}


process TRIMMOMATIC {
    tag "Trimmomaticing"
    // Only published result from proccess is the trimmed fastq.gz's and can be found in ./results/Trimmomatic
    publishDir "${params.publish_dir}/Trimmomatic", mode: 'copy'

    conda "bioconda::trimmomatic=0.39"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
        // important note: These paths do not go to file level but folder unlike other proccesses (for file filtering proccess)
        path shortreads1 // channel of short read pairs 1 PATH/TO/SHORT/READS1
        path shortreads2 // channel of short read pairs 1 PATH/TO/SHORT/READS2
        path truseq // Path to TruSeq3-PE-2.fa required for trimming
    output:
	    path 'trimmomatic'
    
    script:
    """
    mkdir trimmomatic
    cp -r ${shortreads1} trimmomatic
    cp -r ${shortreads2} trimmomatic
    cp ${truseq} trimmomatic
    cd trimmomatic

    trimmomatic PE -threads 40 ${shortreads1}/*_R1_001.fastq.gz ${shortreads1}/*_R2_001.fastq.gz \
        short_reads.1_1.trim.fastq.gz short_reads.1_1un.trim.fastq.gz \
        short_reads.1_2.trim.fastq.gz short_reads.1_2un.trim.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:36 ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15

    trimmomatic PE -threads 40 ${shortreads2}/*_R1_001.fastq.gz ${shortreads2}/*_R2_001.fastq.gz \
        short_reads.2_1.trim.fastq.gz short_reads.2_1un.trim.fastq.gz \
        short_reads.2_2.trim.fastq.gz short_reads.2_2un.trim.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:36 ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15
    """
}