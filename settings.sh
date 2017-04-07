
# fasta file of assembly to scaffold
fafile=FAFILE

# fastq file with long reads
fqfile=FQFILE


#bin directory: where the executables are
bindir=$MYSMISDIR/smissv-bin
# working directory
odir=ODIR

#which samtools
mysamtools=`which samtools`
# which aligner
aligner=bwa


# min length of long reads 
minlen=3000
# scaffold settings
min_edge=5.0
num_sigma=`awk "BEGIN {print ($min_edge)/2.5}"`
# length of fake mate pairs
fakelen=2000
# distance between one fake mate pairs and the next
step=200

### Alignment settings
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


# running mode
debug=0

# Scaffolding settings
# Settings for selecting indels
distance_length=2500
# minimum alignment score 
min_mapscore=10
# minimum identity between long read and reference (%)
minID=60


