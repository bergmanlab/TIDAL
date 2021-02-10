#!/bin/sh

#module the necessary modules of software (grid engine specific), not needed in other OS
#source /projectnb/lau-bumc/SOFTWARE/TIDAL/github/TIDAL1.2/CODE/TIDAL_module.sh
#source /projectnb/lau-bumc/SOFTWARE/TIDAL/CODE/TIDAL_module.sh

#runs the TIDAL pipeline
CODEDIR=$3

#pass the fastq filename as argument
lib=$1
prefix=${1%.fastq*}
#read_len=151
read_len=$2

genomedb=$4
masked_genomedb=$5
consensus_TEdb=$6
GENOME=$7
refseq_annotationfile=$8
chrlen_file=$9
chrDir=${10}
gemMappabilityFile=${11}
fly_virus_structure_repbase_DB=${12}
MASKED_GENOME=${13}
repeat_masker_file=${14}
table_lookup=${15}
TE_fasta=${16}

#data prep and creation of uq file
$CODEDIR/data_prep.sh $lib $CODEDIR 

#run the insertion part of TIDAL
$CODEDIR/insert_pipeline.sh \
	$lib".uq.polyn" \
	$read_len \
	$CODEDIR \
	$genomedb \
	$masked_genomedb \
	$consensus_TEdb \
	$GENOME \
	$refseq_annotationfile \
	$chrlen_file \
	$chrDir \
	$gemMappabilityFile \
	$fly_virus_structure_repbase_DB \
	$TE_fasta

#set up symbolic links to do the depletion part of TIDAL
$CODEDIR/setup.sh $lib
#run the depletion part of TIDAL
$CODEDIR/depletion_pipeline.sh \
	$lib".uq.polyn" \
	$read_len \
	$CODEDIR \
	$genomedb \
	$masked_genomedb \
	$consensus_TEdb \
	$GENOME \
	$MASKED_GENOME \
	$repeat_masker_file \
	$refseq_annotationfile \
	$table_lookup \
	$chrlen_file
	
#compile insertion and depletion results
$CODEDIR/last_part.sh $lib $CODEDIR

#run the TE_coverage_TIDAL.sh 
pushd insertion
$CODEDIR/TE_coverage_TIDAL.sh \
	$lib".uq.polyn" \
	$prefix".sort.bam" \
	$CODEDIR \
	$consensus_TEdb \
	$TE_fasta
popd
