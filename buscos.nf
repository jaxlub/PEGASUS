#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BUSCOS {
    tag "Busco"
    publishDir "${params.publish_dir}/Busco", mode: 'copy' 

    input:
        path final_genome_path // output channel of HAPOGMAP which contains hapog_results/hapog.fasta
        path ragtag_genome // output channel of RAGTAG which contains ragtag_output/ragtag.scaffold.fasta
        path busco // parameter to directory of known busco genome
        path busco_image
    output:
        path 'busco_genome'
        path 'busco_ragtag'
    script:
    """
    module load singularity
    cp ${final_genome_path} \$PWD/results.fasta
    cp ${ragtag_genome}/ragtag.scaffold.fasta \$PWD/ragtag.fasta

    singularity exec -B \$PWD ${busco_image} busco -m genome -i results.fasta -o busco_genome -l ${busco}

    singularity exec -B \$PWD ${busco_image} busco -m genome -i ragtag.fasta -o busco_ragtag -l ${busco}


    """
}

/*
process HAPOG12 {
    tag "Mapping and SAMing to BAMing"

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
        path flye // out channel of flye which is path to flye results containing assembly.fasta
	    path cleaned_short_reads // Out channel of Short Centrifuge which contains folder of .trim.fastq.gz files
    output:
        path 'short_reads_mapped2_sorted.bam' 
    script: 
    """
    # Index FLYE results
    bwa-mem2 index -p mem2 ${flye}

    # Map
    bwa-mem2 mem mem2 -t 40 ${cleaned_short_reads}/clean_short_reads_1.trim.fastq.gz ${cleaned_short_reads}/clean_short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    # Convert to Bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    # Sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam
    """
}

process HAPOG23 {
    tag "HAPOGing"
    publishDir "${params.publish_dir}/Hapog23", mode: 'copy' 

    input:
        path flye // channel output of flye which is assembly.fasta
        path hapog_mapping // channel output of HAPOG1 which is short_reads_mapped2_sorted.bam
        path hapog // path to hapog.py

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    # cd ${hapog_mapping}
    # cd .. # should be in hapog/ now
    # each process is own work folder so might not be worth ^, alternative below
    # mkdir hapog_polish
    # might be able to get away with this below.... work should I have flye,hapog_mapping dir and then make polishing to out channel.
    
    """
    mv hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa ./hapog_result.fasta
    python3 /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/AlphaNumTest/sanitizingHeaders.py hapog_result.fasta
    python3 ${hapog} --genome sanitized3_fasta.fasta -b short_reads_mapped2_sorted.bam -o polishing -t 36 -u
    """
}
*/
params.buscoPath = '/gpfs1/home/e/g/eguswa/scratch/bullhead/actinopterygii_odb10'
params.genome_path = '/users/j/l/jlubkowi/scratch/bullhead_project/manual_run_LB1/ntLink/hapog/polishing/hapog_results/hapog.fasta'

workflow {
BUSCOS(params.genome_path, params.buscoPath)

}