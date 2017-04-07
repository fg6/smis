
Scaffolding pipeline using data from long reads technologies (PacBio, ONT)
to scaffold an initial draft assembly. The long reads are shred in smaller segments 
(f.i. 1000 bp) to create fake mate-pairs. The fake mates are
then aligned against the draft assembly and the spinner scaffolder looks for
links between contigs and create scaffolds. 

Original pipeline from Zemin Ning (zn1@sanger.ac.uk): http://www.sanger.ac.uk/science/tools/smis.
Version here modified to use the bwa aligner instead of smalt, and to automize
compiling and running.

### Download and Compile:
Requirements for compiling: zlib, bamtools

	$ export MYBAMTOOLS=/full/path/to/bamtools
	$ git clone https://github.com/fg6/smis.git
	$ cd smis 
	$ export MYSMISDIR=`pwd`
	$ ./makeall.sh

(Tested with gcc-4.9.2, zlib-1.2.8, bamtools-2.4.0) 

### Run 
#### Setup 
Make sure to have $MYSMISDIR set before launching "setup.sh"

	export MYSMISDIR=`pwd`
	$MYSMISDIR/setup.sh </full/path/to/destdir> <draft_assembly> <long_reads>

	where:
   	   /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
   	   draft_assembly: fasta file of the assembly to be scaffolded
  	   long_reads: fastq file of long reads for scaffolding
   

#### Run:
Requirements for running: samtools, bwa (or smalt) in PATH.

	cd /full/path/to/destdir
   	./mysmissv.sh

(Tested with samtools-1.3.1, bwa-0.7.12)

#### Results

Scaffolds will be in /full/path/to/destdir/spinner_scaffolds.fasta


