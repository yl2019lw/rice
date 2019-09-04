#!/bin/bash
# Author:long wei
# Date:2017-11-12 15:02:29

# This script uses salmon to do rna seq quant

PROJECT=$PWD

#data dir, files can be autodetected
DATA_DIR=/media/longwei/OS/Users/longw/Project/first/data/qylHP 
#TODO:auto detect sample ext
SAMPLE_EXT=(_R1.fastq.gz _R2.fastq.gz) #trim ext to get sample name

#out dir
OUT_DIR=$PROJECT/out

#thread count
THREADS=8

#usage
usage()
{
    echo "Usage: `basename $0` [option]
    Preprocessing HTS data refer to ranseq_hisat.sh which call trimmomatic.

Options:
    data storage option:
        -d  <DIR>   location of HTS data directory, contains fasta/fastq files
        -r  <FILE>  path of reference transcriptome fna/fasta file
        -o  <DIR>   location of output directory, default under project dir

    performance option:
        -n  <INT>   threads count, default 8 

    help option:
        -h  show help
    "
}
#end of usage

#parse argument
parse_opt()
{
    while getopts "d:r:o:n:h" opt
    do
        case "$opt" in
        d)
            DATA_DIR=$OPTARG;;
        r)
            REF_TRANSCRIPTOME=$OPTARG;;
        o)
            OUT_DIR=$OPTARG;;
        n)
            THREADS=$OPTARG;;
        h)
            usage
            exit 0;;
        *)
            usage
            exit 1;;
        esac
    done
}
#end of parse opt

#detect name of samples
detect_samples()
{
    samplef=$OUT_DIR/samples.txt
    ls $DATA_DIR > $samplef
    for ext in ${SAMPLE_EXT[@]}
    do
        sed -i "s/$ext//g" $samplef
    done
    sed -i "s/_unpaired//g" $samplef
    sed -i "/^trim$/d" $samplef
    { rm $samplef && uniq > $samplef; } < $samplef
    samples=`cat $samplef`
}
#end of detect samples

main()
{
    parse_opt $@

    type salmon > /dev/null 2>&1 || { echo "Please install salmon"; exit 1; }

    [ ! -e $REF_TRANSCRIPTOME ] && { echo "Please check ref transcriptome:$REF_TRANSCRIPTOME"; exit 1; }

    [ -z "`find $DATA_DIR -type f`" ] && { echo "Please check data dir:$DATA_DIR"; exit 1; }

    #salmon dir
    SALMON_DIR=$OUT_DIR/salmon
    SALMON_INDEX=$SALMON_DIR/index
    SALMON_QUANT=$SALMON_DIR/quant
    SALMON_INDEX_LOG=$SALMON_DIR/salmon_index.log

    [ ! -d $SALMON_QUANT ] && mkdir -p $SALMON_QUANT

    detect_samples

    salmon index -t $REF_TRANSCRIPTOME -i $SALMON_INDEX > $SALMON_INDEX_LOG 2>&1

    for sample in $samples
    do
        local r1=${DATA_DIR}/$sample${SAMPLE_EXT[0]}
        local r2=${DATA_DIR}/$sample${SAMPLE_EXT[1]}
        salmon quant -i $SALMON_INDEX -l A -1 ${r1} -2 ${r2} -p $THREADS -o $SALMON_QUANT/$sample
    done
}

main $@
