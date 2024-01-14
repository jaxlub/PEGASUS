#!/usr/bin/env nextflow
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


params.truseq3 = "/users/j/l/jlubkowi/scratch/bullhead_project/HL4/TruSeq3-PE-2.fa"
truseq_ch = Channel.fromPath( params.truseq3, checkIfExists: true)

params.centrifugeRscript = "/users/j/l/jlubkowi/scratch/bullhead_project/HL42.0/run_centrifuge_clean.R"
centrifugeRscript = Channel.fromPath( params.centrifugeRscript, checkIfExists: true)

phvDatabase = Channel.fromPath( params.phvDatabase, checkIfExists: true)
shortreads1_ch = Channel.fromPath( params.shortPath1, checkIfExists: true )
shortreads2_ch = Channel.fromPath( params.shortPath2, checkIfExists: true )
longread_ch = Channel.fromPath( params.nanoPath , checkIfExists: true)
quastpathGFF_ch = Channel.fromPath( params.quastPathGFF, checkIfExists: true)
quastpathFNA_ch = Channel.fromPath( params.quastPathFNA, checkIfExists: true)
buscopath_ch = Channel.fromPath( params.buscoPath, checkIfExists: true)
ref1_ch = Channel.fromPath( params.ref1, checkIfExists: true)


process PORECHOP_ABI {
    publishDir "${params.publish_dir}/pore_chop", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop_abi:0.5.0--py310h590eda1_0':
        'biocontainers/porechop_abi:0.5.0--py310h590eda1_0' }"

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        path long_reads

    output:
        path 'long_reads_clean.fastq.gz'

    script:
    """
    porechop_abi --input $long_reads --threads 10 --output long_reads_clean.fastq.gz 
    """
}

process FASTQC {
    publishDir "${params.publish_dir}/FastQC", mode: 'copy'

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
    publishDir "${params.publish_dir}/long_centrifuge", mode: 'copy'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
        path longreads // path for longread
        path database // Path to p+h+v indexes

    output:
        path 'long_cent_out.txt'
        path 'longreads.fastq'

    script:
    """
        pwds=\$(echo \$PWD)
        cp ${longreads} \$pwds/longreads.fastq.gz
        cp ${database}/* \$pwds

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


process SHORT_CENTRIFUGE {
    publishDir "${params.publish_dir}/short_centrifuge", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
        path trimmomatic // path to trimmomatic result
        path database // Path to p+h+v.tar.gz

    output:
        path 'short_cent'
        path 'short_cent_out.txt'
    
    script:
        """
        cat ${trimmomatic}/*.1_1.trim.fastq.gz trimmomatic/*.2_1.trim.fastq.gz > short_reads_1.trim.fastq.gz
        cat	${trimmomatic}/*.1_2.trim.fastq.gz trimmomatic/*.2_2.trim.fastq.gz > short_reads_2.trim.fastq.gz
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
        path short_inputs // path with file containing short_reads
        path short_cent // output .txt of centrifuge of files to remove

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
        path sequences // file short_cent containing shortreads
        path not_contam // file containing not contaniminated reads

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





process NANOPLOT {
    publishDir "${params.publish_dir}/Nanoplot", mode: 'copy'

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

    publishDir "${params.publish_dir}/Flye", mode: 'copy'

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
    // Only published result from proccess is the trimmed fastq.gz's and can be found in ./results/Trimmomatic
    publishDir "${params.publish_dir}/Trimmomatic", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
        // important note: These paths do not go to file level but folder unlike other proccesses (for file filtering proccess)
        path shortreads1 // channel of short read pairs 1 PATH/TO/SHORT/READS1
        path shortreads2 // channel of short read pairs 1 PATH/TO/SHORT/READS2
        path truseq // Path to TruSeq3-PE-2.fa required for trimming
        path porechop_abi // Porechop ABI output to ensure it has resources to run
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

process BUSCOS {
    tag "Busco"
    publishDir "${params.publish_dir}/Busco", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.5.0--pyhdfd78af_0' :
        'biocontainers/busco:5.5.0--pyhdfd78af_0' }"



    input:
        path final_genome_path // output channel of HAPOGMAP which contains hapog_results/hapog.fasta
        path ragtag_genome // output channel of RAGTAG which contains ragtag_output/ragtag.scaffold.fasta
        path busco // parameter to directory of known busco genome
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

process HAPOG {
    publishDir "${params.publish_dir}/Hapog1", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hapog:1.3.7--py39hf95cd2a_0' :
        'biocontainers/hapog:1.3.7--py39hf95cd2a_0' }"

    input:
        path flye // channel output of flye which is assembly.fasta
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    """
    hapog --genome assembly.fasta --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    """
}


process NTLINK {
    tag "ntLinking"
    publishDir "${params.publish_dir}/NTLink", mode: 'copy' 

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ntlink:1.3.9--py39hd65a603_2' :
        'biocontainers/ntlink:1.3.9--py39hd65a603_2' }"

    input:
        path nano_reads // Channel of all long reads PATH/TO/LONG/READ/*
        path hapog_polish // Out channel of HAPOG containing hapog.fasta
    output:
        path 'hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa'
    script: 
    """
    pwds=\$(echo \$PWD)
    cp ${hapog_polish} \$pwds/hapog_result.fasta
    cp ${nano_reads} \$pwds/nanopore_raw.fastq.gz
    ntLink_rounds run_rounds target=hapog_result.fasta reads=nanopore_raw.fastq.gz k=32 w=250 t=16 overlap=True rounds=5
    """
}


process HAPOGMAP {
    tag "HAPOGing the NTLINK"
    publishDir "${params.publish_dir}/Hapog2", mode: 'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hapog:1.3.7--py39hf95cd2a_0' :
        'biocontainers/hapog:1.3.7--py39hf95cd2a_0' }"

    input:
        path ntlink // Out channel of ntlink which contains hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    """
    hapog --genome hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    """
}

process QUASTGENOME {
    tag "Quasting"

    publishDir "${params.publish_dir}/Quast_Genome", mode: 'copy' 

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path final_genome_path // output channel of HAPOGMAP which contains hapog_results/hapog.fasta
        path quastGenomeGFF // parameter of path to file GFF
        path quastGenomeFNA // parameter of path to file 
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

    publishDir "${params.publish_dir}/Quast_Ragtag", mode: 'copy' 

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path ragtag_genome // Output channel of Ragtag which contains ragtag_output/ragtag.scaffold.fasta
        path quastGenomeGFF // parameter of path to file GFF
        path quastGenomeFNA // parameter of path to file 

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
 
process MULTIQC {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"
    
    input:
        path busco
        path mosdepth
        path quast_genome
        path quast_ragtag
        path unknown

    output:
        path 'multiqc_report.html'
        path 'multiqc_data'
    script:
        """
        pwds=\$(echo \$PWD)
        cd ../../.
        multiqc . -d
        cp multiqc_data \$pwds -r
        cp multiqc_report.html \$pwds
        """
}

process BAMCONVERSION {
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
        path genome // Out channel of HAPOGMAP which contains results.fasta
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files
    
    output:
        path 'short_reads_mapped2_sorted.bam' // channel output to bam file
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
        path genome // path to shortreads in bam file
        path index // path to indexes of shortreads in bai file
    output:
        path 'mosdepth_result.*'
    script: 
    
    """
    mosdepth -t 8 mosdepth_result ${genome}
    """
}

process RAGTAG {
    tag "Ragging and Tagging "
    publishDir "${params.publish_dir}/RagTag", mode: 'copy' 
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0' :
        'biocontainers/ragtag:2.1.0--pyhb7b1952_0' }"

    input:
        path genome
        path ref1
    output:
        path 'ragtag_results'
    script: 
    """
    pwds=\$(echo \$PWD)
    cp ${ref1} \$pwds/ref1.fasta
    cp ${genome} \$pwds/results.fasta


    ragtag.py scaffold -o ragtag_results ref1.fasta results.fasta

    """
}

workflow {
    PORECHOP_ABI(longread_ch)


// Stats
NANOPLOT(longread_ch)
FASTQC(shortreads1_ch, shortreads2_ch, longread_ch)

// Short Read prep
TRIMMOMATIC(shortreads1_ch, shortreads2_ch, truseq_ch, PORECHOP_ABI.out)
SHORT_CENTRIFUGE(TRIMMOMATIC.out, phvDatabase)
R_PROCESSING_SHORT(centrifugeRscript, SHORT_CENTRIFUGE.out)
SHORT_SEQTK(R_PROCESSING_SHORT.out)


// Long Read prep
//PORECHOP_ABI(longread_ch)
LONG_CENTRIFUGE(PORECHOP_ABI.out, phvDatabase)
//LONG_CENTRIFUGE(longread_ch, phvDatabase)
R_PROCESSING_LONG(centrifugeRscript, LONG_CENTRIFUGE.out)
LONG_SEQTK(R_PROCESSING_LONG.out)
FLYE(LONG_SEQTK.out)

// first build
HAPOG(FLYE.out, SHORT_SEQTK.out)
NTLINK(longread_ch, HAPOG.out)
HAPOGMAP(NTLINK.out, SHORT_SEQTK.out)
RAGTAG(HAPOGMAP.out, ref1_ch)

// Quality Checks
BUSCOS(HAPOGMAP.out, RAGTAG.out, buscopath_ch)
QUASTGENOME(HAPOGMAP.out,quastpathGFF_ch, quastpathFNA_ch)
QUASTRAGTAG(RAGTAG.out,quastpathGFF_ch, quastpathFNA_ch)
BAMCONVERSION(HAPOGMAP.out, SHORT_SEQTK.out)
MOSDEPTH(BAMCONVERSION.out)

MULTIQC(BUSCOS.out, MOSDEPTH.out, QUASTGENOME.out,QUASTRAGTAG.out)


}

