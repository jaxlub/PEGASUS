#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.publish_dir = './Results'

include { FASTQC } from './submit1process'
include { LONG_CENTRIFUGE } from './submit1process'
include { SHORT_CENTRIFUGE } from './submit1process'
include { NANOPLOT } from './submit1process'
include { FLYE } from './submit1process'
include { TRIMMOMATIC } from './submit1process'


include { HAPOG } from './submit2process'
//include { HAPOG2 } from './submit2process'
include { NTLINK } from './submit2process'
include { HAPOGMAP } from './submit2process'
//include { NTLINKPOLISH } from './submit2process'
//include { NTLINKHAPOG } from './submit2process'


//include { BUSCOFLYE } from './submit3process'
//include { BUSCOHAPOG } from './submit3process'
//include { QUASTFLYE } from './submit3process'
//include { QUASTHAPOG } from './submit3process'
include { QUASTGENOME } from './submit3process'
include { QUASTRAGTAG } from './submit3process'
//include { BUSCO } from './submit3process'
include { MULTIQC } from './submit3process'

include { MOSDEPTH } from './submit4process'
include { RAGTAG } from './submit4process'
include { BAMCONVERSION } from './submit4process'

include { BUSCOS } from './buscos.nf'
//include { HAPOG12 } from './busco_testing/buscos.nf'
//include { HAPOG23 } from './busco_testing/buscos.nf'

params.nanoPath = ""
params.shortPath1 = ""
params.shortPath2 = ""
//params.project_location = ""
params.buscoPath = ""
params.quastPathFNA = ""
params.quastPathGFF = ""
params.trimpath1 = ""
params.trimpath2 = ""
params.hapog = ""

params.truseq3 = "./TruSeq3-PE-2.fa"
truseq_ch = Channel.fromPath( params.truseq3, checkIfExists: true)

params.centrifugeRscript = "./run_centrifuge_clean.R"
centrifugeRscript = Channel.fromPath( params.centrifugeRscript, checkIfExists: true)

params.phvDatabase = "./PHVindexes"
phvDatabase = Channel.fromPath( params.phvDatabase, checkIfExists: true)

params.mosdepth_image = "./containers/AP_singu_mosdepth.sif"
mosdepth_img = Channel.fromPath( params.mosdepth_image, checkIfExists: true)

params.ntlink_image = "./containers/ntlink.sif"
ntlink_img = Channel.fromPath( params.ntlink_image, checkIfExists: true)

params.ragtag_image = "./containers/ragtag.sif"
ragtag_img = Channel.fromPath( params.ragtag_image, checkIfExists: true)

params.r_image = "./containers/rstudio.simg"
r_img = Channel.fromPath( params.r_image, checkIfExists: true)

params.centrifuge_image = "./containers/singularity-recipes_centrifuge.sif"
centrifuge_img = Channel.fromPath( params.centrifuge_image, checkIfExists: true)

params.seqtk_image = "./containers/singularity-seqtk_latest.sif"
seqtk_img = Channel.fromPath( params.seqtk_image, checkIfExists: true)

params.quast_image = "./containers/quast.sif"
//params.quast_image = "./containers/assembly-utils_quast_5.0.2.sif"
quast_img = Channel.fromPath( params.quast_image, checkIfExists: true)

params.busco_image = "./containers/busco.img"
busco_img = Channel.fromPath( params.busco_image, checkIfExists: true)

params.multiqc_image = "./containers/multiqc:1.14--pyhdfd78af_0"
multiqc_img = Channel.fromPath( params.multiqc_image, checkIfExists: true)


params.ref1 = ""

//params.ref2 = ""
//params.ref3 = ""
//params.quastGFF = ""
//params.guastFNA = ""


shortreads1_ch = Channel.fromPath( params.shortPath1, checkIfExists: true )
shortreads2_ch = Channel.fromPath( params.shortPath2, checkIfExists: true )
longread_ch = Channel.fromPath( params.nanoPath , checkIfExists: true)
trimpath1_ch = Channel.fromPath( params.trimpath1, checkIfExists: true)
trimpath2_ch = Channel.fromPath( params.trimpath2, checkIfExists: true)
hapogpath_ch = Channel.fromPath( params.hapog, checkIfExists: true)
quastpathGFF_ch = Channel.fromPath( params.quastPathGFF, checkIfExists: true)
quastpathFNA_ch = Channel.fromPath( params.quastPathFNA, checkIfExists: true)
buscopath_ch = Channel.fromPath( params.buscoPath, checkIfExists: true)
ref1_ch = Channel.fromPath( params.ref1, checkIfExists: true)
//ref2_ch = Channel.fromPath( params.ref2, checkIfExists: true)
//ref3_ch = Channel.fromPath( params.ref3, checkIfExists: true)



workflow {
allshortreads_ch = shortreads1_ch.mix(shortreads2_ch)
// Stats
NANOPLOT(longread_ch)
FASTQC(trimpath1_ch, trimpath2_ch, longread_ch)

// Short Read prep
TRIMMOMATIC(trimpath1_ch, trimpath2_ch, truseq_ch)
SHORT_CENTRIFUGE(centrifugeRscript, TRIMMOMATIC.out, phvDatabase, r_img, centrifuge_img, seqtk_img)

// Long Read prep
LONG_CENTRIFUGE(centrifugeRscript, longread_ch, phvDatabase,r_img, centrifuge_img, seqtk_img)
FLYE(LONG_CENTRIFUGE.out)

// first build
HAPOG(FLYE.out, SHORT_CENTRIFUGE.out, hapogpath_ch)
//HAPOG2(FLYE.out, HAPOG1.out, hapogpath_ch)
NTLINK(longread_ch, HAPOG.out, ntlink_img)
HAPOGMAP(NTLINK.out, SHORT_CENTRIFUGE.out, hapogpath_ch)
RAGTAG(ragtag_img, HAPOGMAP.out, ref1_ch)

// Quality Checks
BUSCOS(HAPOGMAP.out, RAGTAG.out, buscopath_ch, busco_img)
QUASTGENOME(HAPOGMAP.out,quastpathGFF_ch, quastpathFNA_ch, quast_img)
QUASTRAGTAG(RAGTAG.out,quastpathGFF_ch, quastpathFNA_ch, quast_img)
BAMCONVERSION(HAPOGMAP.out, SHORT_CENTRIFUGE.out)
MOSDEPTH(mosdepth_img, BAMCONVERSION.out)

MULTIQC(multiqc_img)

//MULTIQC(multiqc_img, BUSCOS.out, TRIMMOMATIC.out, MOSDEPTH.out, FASTQC.out, QUASTGENOME.out, QUASTRAGTAG.out, NANOPLOT.out)
}

