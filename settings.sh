
# fasta file of assembly to scaffold
fafile=FAFILE

# fastq file with long reads
fqfile=FQFILE

# working directory
odir=ODIR

#bin directory: where the executables are
bindir=$MYSMISDIR/smissv-bin

#which samtools
mysamtools=`which samtools`

# running mode
debug=0



########### FAKE-MATES PARAMETERS: ###########
# min length of long reads 
minlen=3000
# length of fake mates
fakelen=2000
# distance between one fake mate pairs and the next
step=200


########### ALIGNMENT PARAMETERS: ###########

# which aligner to use (bwa, smalt)
aligner=bwa

# bwa
mybwa=`which bwa`
# smalt 
mysmalt=`which smalt`
# max threads:  
nodes=25
#  minimum smith-waterman alignment score to report a hit 
mapscore=50

# Alignment scores and penalties
match=2
subst=-1
gapopen=-1
gapext=-1



########### SCAFFOLDING PARAMETERS: ###########
# minimum numbers of edges:
min_edge=5.0
num_sigma=`awk "BEGIN {print ($min_edge)/2.5}"`

