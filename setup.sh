#!/bin/bash


odir=$1    
fafile=$2
fqfile=$3


if [ $# -lt 3 ]  || [ $1 == '-h' ]; then
        echo; echo "  Usage:" \$MYSMISDIR/$(basename $0) \</full/path/to/destdir\> \<draft_assembly\> \<long_reads\> ; echo
        echo "    /full/path/to/destdir: folder where to run the pipeline (Please provide full path)"
	echo "    draft_assembly: fasta file of the assembly to be scaffolded (Please provide full path)"
	echo "    long_reads: fastq file of long reads for scaffolding (Please provide full path)"
	exit 1
fi

echo; echo "${MYSMISDIR:?Need to set MYSMISDIR to location of your smis}" > /dev/null

mkdir -p $odir
cd $odir

rm -f mysettings.sh
echo export "MYSMISDIR="$MYSMISDIR > mysettings.sh
sed -e "s#ODIR#$odir#g"  $MYSMISDIR/settings.sh \
	| sed -e "s#FAFILE#$fafile#g" | sed -e "s#FQFILE#$fqfile#g"   >> mysettings.sh


sed "s#ODIR#$odir#g" $MYSMISDIR/smissv.sh > mysmissv.sh
chmod +x mysmissv.sh


echo
echo " You're all set! You can now launch $odir/mysmissv.sh to run the pipeline:"
echo "  $ cd "$odir
echo "  $ ./mysmissv.sh"
echo; echo  "  If needed, modify alignment and scaffolding settings in file $odir/mysettings.sh"
echo
