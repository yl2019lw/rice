#!/bin/bash
# Author:long wei
# Date: 2017-11-12 20:56:27
# This script uses hisat + samtools + htseq-count to do rna analysis

PROJECT=$PWD

#data dir, files can be autodetected
DATA_DIR=/media/longwei/OS/Users/longw/Project/first/data/qylHP 
#TODO:auto detect sample ext
SAMPLE_EXT=(_R1.fastq.gz _R2.fastq.gz) #trim ext to get sample name

#ref data, files can be autodetected
#REF_DIR=$PROJECT/ref/rice_annotation
REF_DIR=$PROJECT/ref/ncbi

#out dir
OUT_DIR=$PROJECT/out

#thread count
THREADS=8

#trimmomatic jar & adapter location
TRIMMOMATIC_HOME=$HOME/software/Trimmomatic-0.36

#run step control flag
STEP_ALL="true"
STEP_INDEX="false"
STEP_QUALITY="false"
STEP_PREPROCESS="false"
STEP_MAP="false"
STEP_BAM="false"
STEP_COUNT="false"

#usage
usage()
{
    echo "This script uses hisat + smatools + htseq-count to do RNA-Seq count.

Usage: `basename $0` [option]
    In normal case, only need to specity data storage option.
    Reference dir must contain fasta/fa/fna file, at least one of gff & gtf.
    When need to preprocess HTS data, specify trimmomatic home by -t.
    All steps will be executed by default(include preprocess).
    For debug mode, can specify several steps you want to run.

Options:
    data storage option:
        -d  <DIR>   location of HTS data directory, contains fasta/fastq files
        -r  <DIR>   location of reference data directory, contains ref dna & gff/gtf
        -o  <DIR>   location of output directory, default under project dir

    performance option:
        -n  <INT>   threads count, default 8 

    execute step option:
        -i  run hisat2-build with ref to build index
        -q  run fastqc to check HTS data quality
        -p  run trimmomatic to preprocess origin HTS data
        -m  run hisat2 to map sequence with hisat index(known as align)
        -b  run samtools to convert sam file to sorted bam
        -c  run htseq-count to count reads

    miscellaneous option:
        -t  <DIR>   location of trimmomatic direcotry, contains jar & adapter

    help option:
        -h  show help
    "
}
#end of usage

#parse argument
parse_opt()
{
    while getopts "d:r:o:n:iqpmbct:h" opt
    do
        case "$opt" in
        d)
            DATA_DIR=$OPTARG;;
        r)
            REF_DIR=$OPTARG;;
        o)
            OUT_DIR=$OPTARG;;
        n)
            THREADS=$OPTARG;;
        i)
            STEP_ALL="false"
            STEP_INDEX="true";;
        q)
            STEP_ALL="false"
            STEP_QUALITY="true";;
        p)
            STEP_ALL="false"
            STEP_PREPROCESS="true";;
        m)
            STEP_ALL="false"
            STEP_MAP="true";;
        b)
            STEP_ALL="false"
            STEP_BAM="true";;
        c)
            STEP_ALL="false"
            STEP_COUNT="true";;
        t)
            TRIMMOMATIC_HOME=$OPTARG;;
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

#check exec program
check_exec()
{
    softwares=(hisat2 hisat2-build samtools java fastqc htseq-count)
    for soft in $softwares
    do
        type $soft > /dev/null 2>&1 || { echo "Please install ${soft}"; exit 1; }
    done
}
#end of check exec

#detect ref genomic fasta/gff/gtf files
detect_ref_files()
{
    for f in `find $REF_DIR -maxdepth 1 -type f -name "*.gz"`
    do
        gunzip $f
    done

    REF_GENOME=`find $REF_DIR -maxdepth 1 -type f -regex ".*\.\(fa\|fasta\|fna\|seq\)$" | head -1`

    REF_GFF=`find $REF_DIR -maxdepth 1 -type f -regex ".*\.\(gff\|gff3\)$" | head -1`

    REF_GTF=`find $REF_DIR -maxdepth 1 -type f -name "*.gtf" | head -1`

    if [ X$REF_GTF == X ]; then
        REF_GTF=${REF_GFF/gff/gtf}
        REF_GTF=${REF_GFF/gff3/gtf}
    fi

    echo "REF_GENOME:$REF_GENOME"
    echo "REF_GFF:$REF_GFF"
    echo "REF_GTF:$REF_GTF"
}
#end of detect ref files

#init directory
init_directory()
{
    #fastqc dir
    FASTQC_DIR=$OUT_DIR/fastqc

    #hisat file dir, store hisat index files
    HISAT_DIR=$OUT_DIR/hisat
    HISAT_INDEX=$HISAT_DIR/hisat_index
    HISAT_INDEX_LOG=$HISAT_DIR/hisat_index.log
    SPLICE_SITE_FILE=$HISAT_DIR/genome.ss
    EXON_FILE=$HISAT_DIR/genome.exon
    
    #sam/bam dir, store sample sam/bam files
    SAM_DIR=$OUT_DIR/sam
    
    #count dir, store sample counts
    COUNT_DIR=$OUT_DIR/count

    for dir in $DATA_DIR $REF_DIR
    do
        [ -z "`find $dir -type f`" ] && { echo "Please check data, $dir is empty"; exit 1; }
    done

    for d in $OUT_DIR $FASTQC_DIR $HISAT_DIR $SAM_DIR $COUNT_DIR
    do
        [ ! -d $d ] && mkdir -p $d
    done
}
#end of init directory

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

#prepare running env
init_run_environment()
{
    check_exec
    init_directory
    detect_ref_files
    detect_samples
}
#end of init run environment

#preprocess by trim data
trim_data()
{
    #TODO:add fastqc check
    if [[ $STEP_ALL == "false"  && $STEP_PREPROCESS == "false" ]]; then
        return 0
    fi

    local TRIM_ADAPTER=$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa 
    local TRIM_JAR=$TRIMMOMATIC_HOME/trimmomatic-0.36.jar

    if [[ ! -e $TRIM_JAR  || ! -e $TRIM_ADAPTER ]]; then
        { echo "Please supply $TRIM_JAR & $TRIM_ADAPTER"; return 1; }
    fi

    TRIM_DATA_DIR=$DATA_DIR/trim
    if [ ! -d $TRIM_DATA_DIR ];then
        mkdir -p $TRIM_DATA_DIR
    fi

    local sample=${1}
    #local trimlog=$TRIM_DATA_DIR/${sample}_trim.log
    local i1=$DATA_DIR/${sample}${SAMPLE_EXT[0]} 
    local i2=$DATA_DIR/${sample}${SAMPLE_EXT[1]} 
    local p1=$TRIM_DATA_DIR/${sample}${SAMPLE_EXT[0]} 
    local p2=$TRIM_DATA_DIR/${sample}${SAMPLE_EXT[1]} 
    local up1=$TRIM_DATA_DIR/${sample}_unpaired${SAMPLE_EXT[0]} 
    local up2=$TRIM_DATA_DIR/${sample}_unpaired${SAMPLE_EXT[1]} 
    local params="ILLUMINACLIP:$TRIM_ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    local cmd="java -Xmx8g -Xms4g -jar $TRIM_JAR PE -threads $THREADS -phred33 $i1 $i2 $p1 $up1 $p2 $up2 $params"
    echo $cmd
    eval $cmd
}
#end of trim data

#fastqc quality check
fastqc_quality_check()
{
    if [[ $STEP_ALL == "false"  && $STEP_QUALITY == "false" ]]; then
        return 0
    fi
    #or use awk 'ORS=" "' to join lines
    local seqlist=`find $DATA_DIR -maxdepth 1 -type f -name "*.fastq*" | paste -sd " " -`
    local cmd="fastqc -o $FASTQC_DIR $seqlist"
    echo $cmd 
    eval $cmd
}
#end of fastqc quality check

#hisat2 build index
hisat_build_index()
{
    if [[ $STEP_ALL == "false"  && $STEP_INDEX == "false" ]]; then
        return 0
    fi

    [ -z $REF_GENOME ] && { echo "ref genome file not found"; exit 1; }

    [[ -z $REF_GFF && -z $REF_GTF ]]  && { echo "gff/gtf file not found"; exit 1; }

    if [ ! -e $REF_GTF ]; then
        if [ ! -e $REF_GFF ]; then
            { echo "Please supply gff/gtf file"; return 1; }
        else
            if  type gffread > /dev/null ; then 
                gffread $REF_GFF -T -o $REF_GTF
            else
                { echo "Please install gffread from cufflinks"; exit 1; }
            fi
        fi
    fi
 
    if [ ! -e $SPLICE_SITE_FILE ]; then
        if type extract_splice_sites.py > /dev/null ; then
            extract_splice_sites.py $REF_GTF > $SPLICE_SITE_FILE
        else 
            { echo "Please install extract_splice_sites.py from hisat2"; exit 1; }
        fi
    fi

    if [ ! -e $EXON_FILE ]; then
        if type extract_exons.py  > /dev/null ; then
            extract_exons.py $REF_GTF > $EXON_FILE
        else
           { echo "Please install extract_exons.py from hisat2"; exit 1; }
        fi
    fi

    local cmd="hisat2-build -p $THREADS --ss $SPLICE_SITE_FILE --exon $EXON_FILE $REF_GENOME $HISAT_INDEX > $HISAT_INDEX_LOG 2>&1"
    echo $cmd
    eval $cmd
}
#end of hisat build index

#hisat align
hisat_align()
{
    if [[ $STEP_ALL == "false"  && $STEP_MAP == "false" ]]; then
        return 0
    fi

    local sample=$1
    if [ -z "$TRIM_DATA_DIR" ]; then
        ADATA_DIR=$DATA_DIR
    else
        ADATA_DIR=$TRIM_DATA_DIR 
    fi
    local r1=${ADATA_DIR}/$sample${SAMPLE_EXT[0]}
    local r2=${ADATA_DIR}/$sample${SAMPLE_EXT[1]}
    local sam=${SAM_DIR}/$sample.sam
    local log=${SAM_DIR}/$sample.log
    local cmd="hisat2 -p $THREADS -q -x $HISAT_INDEX -1 ${r1} -2 ${r2} -S ${sam} > ${log} 2>&1"
    echo $cmd
    eval $cmd
}
#hisat align end

#convert sam to sorted bam
sam_to_sorted_bam()
{
    if [[ $STEP_ALL == "false"  && $STEP_BAM == "false" ]]; then
        return 0
    fi

    local sample=$1
    local samf=${SAM_DIR}/$sample.sam
    local bamf=${samf/\.sam/\.unsort.bam}
    local sortf=${samf/\.sam/\.bam}

    echo "convert sam file:$samf to sorted bam:$sortf"

    samtools view -@ $THREADS -b -S -o $bamf $samf
    rm $samf > /dev/null
    samtools sort -n -@ $THREADS -O BAM -o $sortf $bamf 
    rm $bamf > /dev/null
}
#end of sam to sorted bam

#htseq-count
htseq_count()
{
    if [[ $STEP_ALL == "false"  && $STEP_COUNT == "false" ]]; then
        return 0
    fi

    local sample=$1
    local bamf=${SAM_DIR}/$sample.bam
    local countf=${COUNT_DIR}/$sample.count
    local count_log=${COUNT_DIR}/$sample.log

    #default id attr is gene_id, for NCBI rice use gene
    local cmd="samtools view -@ $THREADS -h ${bamf} | htseq-count -s no -t gene -i ID - ${REF_GFF} > ${countf} 2>${count_log}"

    echo ${cmd}
    eval ${cmd}
}
#end of htseq count

main()
{
    parse_opt $@
    init_run_environment

    #step fastqc quality check
    fastqc_quality_check

    #step index
    hisat_build_index

    for sample in $samples
    do
        #step preprocess
        trim_data $sample

        #step map
        hisat_align $sample

        #step bam
        sam_to_sorted_bam $sample

        #step count
        htseq_count $sample
    done
}

main $@
