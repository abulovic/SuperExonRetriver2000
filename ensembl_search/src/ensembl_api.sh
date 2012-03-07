#!/bin/bash

ENSEMBLDB="/home/intern/mnt/release-66/fasta/"

if [ $1 == "--help" ]; then
	echo -e "\033[1m\nNAME\n\033[0m"
	echo -e "\t\033[1mensembl_api\033[0m - a command line utility to query a local Ensembl database" 

	echo -e "\033[1m\nSYNOPSIS\n\033[0m"
	echo -e "\tensembl_api --species \033[4mensembl_compatible_species_name\033[0m --id \033[4munique_id\033[0m --id_type \033[4mID_type\033[0m [--strand] [--masked] [--seq_begin] [--seq_end]" | fmt 
	
	echo -e "\033[1m\nDESCRIPTION\n\033[0m"
	echo -e "\tQuery the local Ensembl database for genomic data.\n"
# The mandatory parameters:
	echo -e "\tMandatory arguments are:\n"

	echo -e "\t\033[1m-s --species\033[0m"
	echo -e "\t\tensembl compatible species name (in Latin)"

	echo -e "\t\033[1m-t --id_type\033[0m"
	echo -e "\t\tSequence ID type as specified in the Ensembl FASTA documentation (supercontig, contig, chromosome, clone, scaffold, chunk)" | fmt

	echo -e "\t\033[1m-i --usid\033[0m"
	echo -e "\t\tUnique sequence ID provided by Ensembl"

# Non-mandatory parameters:
	echo -e "\n\tOther possible arguments are:\n"

	echo -e "\t\033[1m-o --output\033[0m"
	echo -e "\t\tThe output file. If it isn't specified, the result will be written in the species_name.fa in the current directory" | fmt
	
	echo -e "\t\033[1m-c --cache\033[0m"
	echo -e "\t\tCached data folder. Contains appropriate seqlevel data. If specified, the search for seqlevel sequences will be performed there." | fmt

	echo -e "\t\033[1m-b --seq_begin\033[0m"
	echo -e "\t\tThe starting point of the desired part of the sequence" | fmt

	echo -e "\t\033[1m-e --seq_end\033[0m"
	echo -e "\t\tThe ending point of the desired part of the sequence" | fmt

	echo -e "\t\033[1m-r --strand\033[0m"
	echo -e "\t\tForward strand: 1 (default), backward strand: -1" | fmt

	echo -e "\t\033[1m-u --unmasked\033[0m"
	echo -e "\t\tUnmasked DNA used if specified, masked default" | fmt

	echo -e "\t\033[1m-a --assembly_version\033[0m"
	echo -e "\t\tThe version of assembly. Used to makes the search faster" | fmt

	echo -e "\033[1m\nINFO\n\033[0m"
	echo -e "\tSoftware freely available.\n\tAuthor: Ana Bulovic\n\n"
	
	exit 0

fi

SPECIES_NAME="";
USID="";
ID_TYPE="";
SEQ_BEGIN=0;
SEQ_END=0;
STRAND=1;
MASKED="T";
ASSEMBLY=""
OUTPUT=""
CACHE_FOLDER=""

# parse the command line options

set -- $(getopt s:t:i:b:e:r:ua:o:c: "$@")

while [ $# -gt 0 ]
do
    case "$1" in

    (-s) SPECIES_NAME=$2; shift;;

    (-t) ID_TYPE=$2; shift;;

    (-i) USID=$2; shift;;

    (-b) SEQ_BEGIN=$2; shift;;

	(-e) SEQ_END=$2; shift;;

	(-r) STRAND=$2; shift;;

	(-u) MASKED=$2; shift;;

	(-a) ASSEMBLY=$2; shift;;

	(-o) OUTPUT=$2; shift;;
	
	(-c) CACHE_FOLDER=$2; shift;;

	(--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*)  break;;
    esac
    shift
done

# check wether any of mandatory parameters are missing

if [[ "$SPECIES_NAME" == "" || "$USID" == "" || "$ID_TYPE" == "" ]] ; then
	echo -e "\033[1m\nERROR\n\033[0m"
	echo -e "\tOne of the mandatory arguments was not specified. Try reading the help by invoking \033[1mensembl_api --help\033[0m.\n\n" | fmt
	exit -1
fi

# if there is no output file specified, generate output file
if [ -z $OUTPUT ]; then
	OUTPUT=$SPECIES_NAME".fa"
fi

# if the sequence BEGINNING and ending have been specified, 
# convert them to modulo/division friendly format ( - 1 )
if [ $SEQ_BEGIN -ne 0 ] && [ $SEQ_END -ne 0 ]; then
	SEQ_BEGIN=$(($SEQ_BEGIN-1))
	SEQ_END=$(($SEQ_END-1))
	echo "Transforming sequence..."
fi

# generate file name from which to search

# if there's a cached folder specified, and not a chromosome in question
if [ -n "$CACHE_FOLDER" ] && [ "$ID_TYPE" != "chromosome" ]; then
	FILE_NAME="$CACHE_FOLDER/$SPECIES_NAME.$ID_TYPE.fa"
else

	#concatenate dna to the folder of the species (lowercase - hence the awk)
	DIR="$ENSEMBLDB/`echo $SPECIES_NAME | awk '{print tolower($0)}'`/dna/"

	#create basic file name - regardless of id type recquired
	#avoid first file - readme
	#TODO : fix this, make more robust
	FILE_NAME=$DIR`ls -1 $DIR | tail -n -2 | head -n 1 | sed -r 's/(.*)dna.*/\1dna/'`

	#with regard to the sequence type, create the file name
	if [ $MASKED == "T" ]; then
		FILE_NAME=$FILE_NAME"_rm"
	fi

	#  with regard to the id type, concatenate chromosome or toplevel
	if [ $ID_TYPE == "chromosome" ]; then
		FILE_NAME=$FILE_NAME"."$ID_TYPE"."$USID".fa"
	else
		FILE_NAME=$FILE_NAME".toplevel.fa"
	fi
	
fi



echo "Reading file "$FILE_NAME"...\n"

# if dealing with chromosome, no search is necessary.

if [ $ID_TYPE == "chromosome" ]; then

	echo "Type of search:	CHROMOSOME"
	# if there are no BEGINNING / ending specified, just copy the file
	if [ $SEQ_BEGIN -eq 0 ] && [ $SEQ_END -eq 0 ]; then
		echo "Length:		WHOLE SEQ"

		echo ">$USID $ID_TYPE:$ASSEMBLY:$USID:whole:$STRAND" > $OUTPUT
		
		# check wether it's the strand or reverse strand
		if [ $STRAND -eq 1 ]; then
			tail -n +2 $FILE_NAME >> $OUTPUT
		else
			tail -n +2 $FILE_NAME | tr -d "\n" | tr [ATGC] [TACG] | rev | fold -b60 >> $OUTPUT
		fi
		echo "output::: chromosome whole"

	# otherwise, extract the data
	else
	
		# check if the coordinates are within boundaries
		
		NAME=`head -n 1 $FILE_NAME`
		SEQ_LEN=`echo $NAME | sed 's/.*:.*:.*:\(.*\):.*/\1/'`
		if [ $SEQ_BEGIN -lt 1 ]; then
			echo -e "Bad coordinates for $ID_TYPE $USID.\nLeft coordinate will be scaled to 1."
			SEQ_BEGIN=1
		fi
		if [ $SEQ_END -gt $SEQ_LEN ]; then
			echo -e "Bad coordinates for $ID_TYPE $USID.\nRight coordinate will be scaled to $SEQ_LEN."
			SEQ_END=$SEQ_LEN
		fi

		BEGINNING=$(($SEQ_BEGIN/60 + 1))	# corresponds to line in the file  where the
											# sequence starts (+1 marks the sequence name line)
		ENDING=$(($SEQ_END/60 + 1))			# here the sequence ends
		echo "Border lines:  ("$BEGINNING")   ("$ENDING")"

		# extract all the necessary lines as a single line
		LINE=`tail -n "+"$(($BEGINNING+1)) $FILE_NAME | head -n $(($ENDING-$BEGINNING+1)) | tr -d "\n"`
	
		# write the whole thing to a file"
		echo ">$USID $ID_TYPE:$ASSEMBLY:$USID:$(($SEQ_BEGIN+1)):$(($SEQ_END+1)):$STRAND" > $OUTPUT

		# get one-line sequence
		LINE=`echo $LINE | tail -c "+"$(($SEQ_BEGIN%60+1)) | head -c $(($SEQ_END - $SEQ_BEGIN + 1))`
		
		# check wether it's the strand or reverse strand
		if [ $STRAND -eq 1 ]; then
			echo $LINE | fold -b60 >> $OUTPUT
		else
			# if reversed, complement and reverse, fold and write to output
			echo $LINE | tr [ATGC] [TACG] | rev | fold -b60 >> $OUTPUT
		fi
		echo "output::: chromosome part"
	fi

# otherwise, if we are not dealing with a chromosome, 
# we have to dig throught toplevel file
else 
	
	echo "Type of search:	$ID_TYPE"	

	# generate regex and find the first line of the entry
	GREPEX=">$USID.*:$ID_TYPE $ID_TYPE:$ASSEMBLY:$USID:.*"
	GREPRES1=`grep -n -s -m 1 -h -E $GREPEX $FILE_NAME`
	# find the line number of the requested entry
	LINE_NUM1=`echo $GREPRES1 | sed -r 's/([0-9]*).*/\1/'`

	#find the next sequence entry (to calculate the end of the sequence)
	GREPRES2=`tail -n "+$(($LINE_NUM1+1))" $FILE_NAME | grep -n -m 1 ">.*"`	# we are searching for the first next occurance of ">"

	# find the line number of next sequence
	LN=`echo $GREPRES2 | sed -r 's/([0-9]*).*/\1/'`	# relative line number
	#LINE_NUM2=$(($LINE_NUM1 + $LN))				# absolute line number
	
	# if there are no beginning / ending specified, just copy the whole sequence to a file
	echo "$SEQ_BEGIN ... $SEQ_END"
	if [ $SEQ_BEGIN -eq 0 ] && [ $SEQ_END -eq 0 ]; then
		echo "Length:		WHOLE SEQ"
		
		# reverse strand?
		if [ $STRAND -eq 1 ]; then
			tail -n "+$LINE_NUM1" $FILE_NAME | head -n $LN > $OUTPUT
		else
			NAME=`tail -n "+$LINE_NUM1" $FILE_NAME | head -n 1`
			NAME=`echo $NAME | sed 's/\(.*\):1/\1:-1'`
			echo $NAME > $OUTPUT
			tail -n "+$(($LINE_NUM1+1))" $FILE_NAME | head -n $(($LN-1)) | tr -d "\n" | tr [ATGC] [TACG] | rev | fold -b60 >> $OUTPUT
		fi
		echo "output::: other whole"

	else
		
		# check if the coordinates are within boundaries
		SEQ_LEN=`echo $GREPRES1 | sed 's/.*:.*:.*:\(.*\):.*/\1/'`
		if [ $(($SEQ_BEGIN+1)) -lt 1 ]; then
			echo -e "Bad coordinates for $ID_TYPE $USID.\nLeft coordinate will be scaled to 1."
			SEQ_BEGIN=1
		fi
		if [ $(($SEQ_END+1)) -gt $SEQ_LEN ]; then
			echo -e "Bad coordinates for $ID_TYPE $USID.\nRight coordinate will be scaled to $SEQ_LEN."
			SEQ_END=$SEQ_LEN
		fi
		
		# write the fasta id to file
		echo "$USID $ID_TYPE:$ASSEMBLY:$USID:$(($SEQ_BEGIN+1)):$(($SEQ_END+1)):$STRAND" > $OUTPUT
		echo "Length:		$SEQ_BEGIN - $SEQ_END"

		# calculate beginning in the sequence
		BEGINNING=$(($SEQ_BEGIN/60+1))

		# calculate the ending
		ENDING=$(($SEQ_END/60+1))
		
		# extract the needed lines
		LINE=`tail -n "+"$(($BEGINNING+$LINE_NUM1)) $FILE_NAME | head -n $(($ENDING-$BEGINNING+1)) | tr -d "\n"`

		# write name to file
		echo ">$USID $ID_TYPE:$ASSEMBLY:$USID:$(($SEQ_BEGIN+1)):$(($SEQ_END+1)):$STRAND" > $OUTPUT
		# remove the uneccessary characters from beginning / end
		LINE=`echo $LINE | tail -c "+"$(($SEQ_BEGIN%60+1)) | head -c $(($SEQ_END-$SEQ_BEGIN + 1))`
		
		# reverse strand?
		if [ $STRAND -eq 1 ]; then
			echo $LINE | fold -b60 >> $OUTPUT
		else
			echo $LINE | tr [ATGC] [TACG] | rev | fold -b60 >> $OUTPUT
		fi
		echo "output::: other part"
	fi
	
fi