#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BAMCONVERSION {
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
        path genome // Out channel of NTLINKHAPOG which contains results.fasta
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files
    
    output:
        path 'short_reads_mapped2_sorted.bam' // channel output to polishing directory containing all results.

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
    """
}
process MOSDEPTH {
    tag "Mosdepthing"
    publishDir "${params.publish_dir}/Mosdepth", mode: 'copy' 

    input:
	    path mosdepth_image // path to the mosdepth image
        path genome // path to genome from second round of hapog
    output:
        path 'mosdepth_result.*'
    script: 
    /*
    Notes: check where and what the samtools file for indexing is coming from
    samtools index ${ntlink_polish}/short_reads_mappedntLink_allMale_genome_sorted.bam 

    # Index FLYE results
    bwa-mem2 index -p mem2 assembly.fasta

    # Map
    bwa-mem2 mem mem2 -t 40 ${cleaned_short_reads}/clean_short_reads_1.trim.fastq.gz ${cleaned_short_reads}/clean_short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    # Convert to Bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    # Sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam

    samtools index ragtag_input.bam
    cp ${ragtag_genome} \$PWD/ragtag_input.bam 
    singularity exec --bind \$PWD ${mosdepth_image} mosdepth -t 8 mosdepth_results ragtag_input.bam 
    */
    """
    module load singularity
    cp ${genome} \$PWD/genome_input.bam 
    samtools index genome_input.bam 
    singularity exec --bind \$PWD ${mosdepth_image} mosdepth -t 8 mosdepth_result genome_input.bam
    """
}

process RAGTAG {
    tag "Ragging and Tagging "
    publishDir "${params.publish_dir}/RagTag", mode: 'copy' 

    input:
	    path ragtag_image // path to the ragtag image
        path genome
        path ref1
    output:
        path 'ragtag_results'
    script: 
    """
    module load singularity
    pwds=\$(echo \$PWD)
    cp ${ref1} \$pwds/ref1.fna
    cp ${genome} \$pwds/results.fasta


    singularity exec --bind \$pwds ${ragtag_image} ragtag.py scaffold -o ragtag_results ref1.fna results.fasta

    """
    //singularity exec --bind \$pwds ${ragtag_image} ragtag.py merge results.fasta out_*/*.agp ragtag_result.map.agp

    /*
input:
	    path ragtag_image // path to the ragtag image
        path genome
        path ref1
        path ref2
        path ref3
    output:
        path 'ragtag_result.map.agp'
    script: 
    """
    pwds=\$(echo \$PWD)
    cp ${ref1} \$pwds/ref1.fna
    cp ${ref2} \$pwds/ref2.fna
    cp ${ref3} \$pwds/ref3.fna


    singularity exec --bind \$pwds ${ragtag_image} ragtag.py scaffold -o out_1 ref1.fna ${genome}
    singularity exec --bind \$pwds ${ragtag_image} ragtag.py scaffold -o out_2 ref2.fna ${genome}
    singularity exec --bind \$pwds ${ragtag_image} ragtag.py scaffold -o out_3 ref3.fna ${genome}

    singularity exec --bind \$pwds ${ragtag_image} ragtag.py merge ${genome} out_*\*.agp ragtag_result.map.agp   
     */
}
