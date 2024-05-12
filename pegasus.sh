#!/bin/bash

nanoPath=""
shortPath1=""
shortPath2=""
quastPathGFF=""
quastPathFNA=""
buscoPath=""
ref1=""
rpath=""
phv=""

usage() {
    echo "Usage: $0 [-n <value>] [-s1 <value>] [-s2 <value>] [-qg <value>] [-qf <value>] [-b <value>] [-r <value>] [-phv <value>]"
    echo "Options:"
    echo "  -n <value>: Path to Nanopore Long Reads"
    echo "  -s1 <value>: Path to Short Reads 1"
    echo "  -s2 <value>: Path to Short Reads 2"
    echo "  -qg <value>: Path to Quast GFF file"
    echo "  -qf <value>: Path to Quast FNA file"
    echo "  -b <value>: Path to BUSCO Reference"
    echo "  -r <value>: Path to reference genome for Ragtag"
    echo "  -phv <value>: Path to PHV files for Centrifuge"
    exit 1
}

check_command() {
    if command -v "$1" &>/dev/null; then
        echo "$1 is installed."
    else
        echo "$1 is not installed."
        echo "Please download and install $1 to proceed."
        exit 1
    fi
}

check_singularity() {
    echo "Checking for Singularity..."
    check_command "singularity"
}

check_nextflow() {
    echo "Checking for Nextflow..."
    check_command "nextflow"
}

check_singularity
check_nextflow
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -n)
            nanoPath="$2"
            shift
            shift
            ;;
        -s1)
            shortPath1="$2"
            shift
            shift
            ;;
        -s2)
            shortPath2="$2"
            shift
            shift
            ;;
        -qg)
            quastPathGFF="$2"
            shift
            shift
            ;;
        -qf)
            quastPathFNA="$2"
            shift
            shift
            ;;
        -b)
            buscoPath="$2"
            shift
            shift
            ;;
        -r)
            ref1="$2"
            shift
            shift
            ;;  
         -phv)
            phv="$2"
            shift
            shift
            ;;  
        *)
            usage
            ;;
    esac
done

script_dir=$(dirname "$0")

pegasus_nextflow=$(echo "$script_dir/bin/pegasus.nf")
rpath=$(echo "$script_dir/bin/run_centrifuge_clean.R")

if [[ -n $nanoPath && -n $shortPath1 && -n $shortPath2 && -n $quastPathGFF && -n $quastPathFNA && -n $buscoPath && -n $ref1 && -n $phv ]]; then
    echo "Configuration 1: Long- and Short-Read Assembly"
    parameters=$(echo "-resume 
    --nanoPath $nanoPath
    --shortPath1 $shortPath1
    --shortPath2 $shortPath2
    --quastPathGFF $quastPathGFF
    --quastPathFNA $quastPathFNA
    --buscoPath $buscoPath
    --ref1 $ref1
    --centrifugeRscript $rpath
    --phvDatabase $phv")

    execute=$(echo "nextflow run $pegasus_nextflow $parameters")
    $execute
else
    echo "Illegal configuration: Incomplete flags entered"
    usage
fi






