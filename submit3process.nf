#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* THESE PROCCESS ARE USED TO CHECK INTERMEDIATE RESULT QAULITY 
process BUSCOFLYE {
    tag "Busco Flye only"

    conda "bioconda::busco=5.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.3--pyhdfd78af_0':
        'biocontainers/busco:5.4.3--pyhdfd78af_0' }"

    input:
        path project_location
        path busco
    script:
    """
    mkdir ${project_location}/busco_flye
    cd ${project_location}/busco_flye

    busco -m genome -i ${project_location}/flye/assembly.fasta -o ${project_location}/busco_flye -l ${busco}
    """
}

process BUSCOHAPOG {
    tag "Busco Hapog only"
    
    conda "bioconda::busco=5.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.3--pyhdfd78af_0':
        'biocontainers/busco:5.4.3--pyhdfd78af_0' }"

    input:
        path project_location
        path busco
    script:
    """
    mkdir ${project_location}/busco_hapog
    cd ${project_location}/busco_hapog

    busco -m genome -i ${project_location}/hapog/polishing/hapog_results/hapog.fasta -o ${project_location}/busco_hapog -l ${busco}
    """
}

process QUASTFLYE {
    tag "Quast Flye"

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path project_location
        path quastGenome
    script:
    """
    mkdir ${project_location}/quast_flye
    cd ${project_location}/quast_flye

    quast ${project_location}/flye/assembly.fasta  -r ${quastGenome}/*.fna -g ${quastGenome}/*.gff
    """
}

process QUASTHAPOG {
    tag "Quast Hapog"

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
        path project_location
        path quastGenome
    script:
    """
    mkdir ${project_location}/quast_hapog
    cd ${project_location}/quast_hapog

    quast ${project_location}/hapog/polishing/hapog_results/hapog.fasta  -r ${quastGenome}/*.fna -g ${quastGenome}/*.gff
    """
}
*/
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
        path quast_image
    output: 
        path 'quast_results'
    script:
    """
    module load singularity
    cp ${quastGenomeGFF} \$PWD/input.gff
    cp ${quastGenomeFNA} \$PWD/input.fna
    cp ${final_genome_path} \$PWD/genome.fasta
    
    singularity exec --bind \$PWD ${quast_image} quast.py genome.fasta -r input.fna -g input.gff --large -t 10
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
        path quast_image
    output: 
        path 'ragtag_quast_results'
    script:
    """
    module load singularity
    cp ${quastGenomeGFF} \$PWD/input.gff
    cp ${quastGenomeFNA} \$PWD/input.fna
    cp ${ragtag_genome}/ragtag.scaffold.fasta \$PWD/ragtag.fasta
        
    singularity exec --bind \$PWD ${quast_image} quast.py ragtag.fasta -r input.fna -g input.gff --large -t 10

    cp -r quast_results ragtag_quast_results
    """
}

/*
process BUSCO {
    tag "Busco"

    publishDir "${params.publish_dir}/Busco", mode: 'copy' 

    input:
        path final_genome_path // output channel of ntlinkpolish which contains hapog_results/hapog.fasta
        path busco // parameter to directory of known busco genome
    output:
        path 'busco'
        
    script:
    """
    module load singularity
    pwds=\$(echo \$PWD)
    cp ${final_genome_path} \$pwds/results.fasta

    singularity exec --bind \$pwds /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/Busco2.sif busco -m genome -i results.fasta -o busco -l ${busco}
    """  
}
*/
/*
    singularity exec bash -c "source activate test_env && pbsv --version"

    singularity exec --bind \$pwds:/mnt /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/Busco2.sif /mnt/busco.sh

    singularity exec --bind \$pwds /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/Busco2.sif busco.sh

    singularity exec --bind \$pwds /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/Busco2.sif "if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
        fi && busco -m genome -i results.fasta -o busco -l ${busco}"

singularity exec --bind \$pwds /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/Busco2.sif bash -c "if [ ! -w \"\$AUGUSTUS_CONFIG_PATH\" ]; then
    # Create writable tmp directory for augustus
    AUG_CONF_DIR=\$( mktemp -d -p /mnt )
    cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
    export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
    echo \"New AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH\"
    fi && busco -m genome -i /mnt/results.fasta -o /mnt/busco -l ${busco}"

*/
// Unsure how to multiqc everything 
process MULTIQC {

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"
    
    input:

        path multiqc_image
        /*
        path busco
        path trimmomatic
        path mosdepth
        path fastqc
        path quast_genome
        path quast_ragtag
        path nanoplot
        path unknown
    */
    output:
        path 'multiqc_report.html'
    script:
        """
        module load singularity
        pwds=\$(echo \$PWD)
        cd ../../.
        singularity exec --bind \$PWD \$pwds/${multiqc_image} multiqc . -d
        cp multiqc_data \$pwds -r
        cp multiqc_report.html \$pwds
        """
}
