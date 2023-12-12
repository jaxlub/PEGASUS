#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process HAPOG {
    publishDir "${params.publish_dir}/Hapog1", mode: 'copy' 

    input:
        path flye // channel output of flye which is assembly.fasta
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files
        path hapog // path to hapog.py

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    """
    python3 ${hapog} --genome assembly.fasta --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    """
}


/*
process HAPOG1 {
    tag "Mapping and SAMing to BAMing"

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
	    path cleaned_short_reads // Out channel of Short Centrifuge which contains folder of .trim.fastq.gz files
        path flye // out channel of flye which is path to flye results containing assembly.fasta
    output:
        path 'short_reads_mapped2_sorted.bam' 
    script: 
    """
    # Index FLYE results
    bwa-mem2 index -p mem2 assembly.fasta

    # Map
    bwa-mem2 mem mem2 -t 40 ${cleaned_short_reads}/clean_short_reads_1.trim.fastq.gz ${cleaned_short_reads}/clean_short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    # Convert to Bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    # Sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam
    """
}
*/
/*
mkdir -p hapog_mapping
    mv trimmomatic hapog_mapping
    cd hapog_mapping


    cat trimmomatic/*.1_1.trim.fastq.gz trimmomatic/*.2_1.trim.fastq.gz > short_reads_1.trim.fastq.gz
    cat	trimmomatic/*.1_2.trim.fastq.gz trimmomatic/*.2_2.trim.fastq.gz > short_reads_2.trim.fastq.gz
 

    # Index FLYE results
    bwa-mem2 index -p mem2 assembly.fasta


    # Map
    bwa-mem2 mem mem2 -t 40 short_reads_1.trim.fastq.gz short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    # Convert to Bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    # Sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam
    """
*/
/*
process HAPOG2 {
    tag "HAPOGing"
    publishDir "${params.publish_dir}/Hapog2", mode: 'copy' 

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
    python3 ${hapog} --genome assembly.fasta -b short_reads_mapped2_sorted.bam -o polishing -t 36 -u
    """
}
*/

process NTLINK {
    tag "ntLinking"
    publishDir "${params.publish_dir}/NTLink", mode: 'copy' 

    input:
        path nano_reads // Channel of all long reads PATH/TO/LONG/READ/*
        path hapog_polish // Out channel of HAPOG2 containing hapog.fasta
        path ntlinkImage // path to image for ntlink env
    output:
        path 'hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa'
    script: 
    """
    module load singularity
    pwds=\$(echo \$PWD)
    cp ${hapog_polish} \$pwds/hapog_result.fasta
    cp ${nano_reads} \$pwds/nanopore_raw.fastq.gz
    singularity exec -B \$pwds \
        ${ntlinkImage} \
        ntLink_rounds run_rounds target=hapog_result.fasta reads=nanopore_raw.fastq.gz k=32 w=250 t=16 overlap=True rounds=5
    
    """
}
//   ln -s ${project_location}/hapog/polishing/hapog_results/hapog.fasta hapog.fasta

/*
process NTLINKPOLISH {
    tag "Polishing the ntLink"
    publishDir "${params.publish_dir}/FinalGenome", mode: 'copy' 


    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
        path ntlink // Out channel of ntlink which contains hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files
    
    output:
        path 'short_reads_mapped2_sorted.bam' // channel output to polishing directory containing all results.

    script: 
    """
    # index
    bwa-mem2 index -p mem2 ${ntlink}

    # map
    bwa-mem2 mem mem2 -t 40 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz > short_reads_mapped2.sam

    #convert to bam
    samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam

    #sort
    samtools sort short_reads_mapped2.bam > short_reads_mapped2_sorted.bam
    """
}
*/
/*
# index
#bwa-mem2 index -p mem2 /users/j/l/jlubkowi/scratch/bullhead_project/manual_run_LB1/ntLink/hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa

# map
#mkdir mappings
#bwa-mem2 mem mem2 LB1_1.trim.fastq.gz LB1_2.trim.fastq.gz > short_reads_mapped2.sam

#comvert to bam
#samtools view -bS short_reads_mapped2.sam > short_reads_mapped2.bam
#rm /users/e/g/eguswa/scratch/bullhead/LB1/hapog/mapping/short_reads_mapped2.sam

#sort
#samtools sort short_reads_mapped2.bam >  short_reads_mapped2_sorted
#rm /users/e/g/eguswa/scratch/bullhead/hapog/SP2_polish-1/mappings/short_reads_mapped2.bam

cd ../

python3 /users/j/l/jlubkowi/scratch/bullhead_project/hapog/hapog.py --genome /users/j/l/jlubkowi/scratch/bullhead_project/manual_run_LB1/ntLink/hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa -b /users/j/l/jlubkowi/scratch/bullhead_project/manual_run_LB1/ntLink/hapog/mapping/short_reads_mapped2_sorted.bam -o polishing -t 16 -u

*/

/*
process NTLINKHAPOG {
    tag "HAPOGing the NTLINK"
    publishDir "${params.publish_dir}/Hapogntlink", mode: 'copy' 

    input:
        path ntlink // Out channel of ntlink which contains hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa
        path polish_out // Out channel of ntlink polish containing short_reads_mapped2_sorted.bam
        path hapog // path to hapog.py

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    """
    mv hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa ./hapog_result.fasta
    python3 ${hapog} --genome /users/j/l/jlubkowi/scratch/bullhead_project/pipeline/AlphaNumTest/sanitized2_fasta.fasta -b short_reads_mapped2_sorted.bam -o polishing -t 16 -u
    """
}
*/


process HAPOGMAP {
    tag "HAPOGing the NTLINK"
    publishDir "${params.publish_dir}/Hapog2", mode: 'copy' 

    input:
        path ntlink // Out channel of ntlink which contains hapog.fasta.k32.w250.z1000.ntLink.5rounds.fa
        path centrifuge_out // Out channel of Short centrifuge which contains hapog/mapping/X_X.trim.fastq.gz files
        path hapog // path to hapog.py

    output:
        path 'polishing/hapog_results/hapog.fasta' // channel output to polishing directory containing all results.

    script: 
    """
    python3 ${hapog} --genome hapog_result.fasta.k32.w250.z1000.ntLink.5rounds.fa --pe1 ${centrifuge_out}/clean_short_reads_1.trim.fastq.gz --pe2 ${centrifuge_out}/clean_short_reads_2.trim.fastq.gz -o polishing -t 16 -u
    """
}