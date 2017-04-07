
Scaffolding pipeline using data from long reads technologies (PacBio, ONT)
to scaffold an initial draft assembly. The long reads are shred in smaller segments 
(f.i. 1000 bp) to create fake mate-pairs. The fake mates are
then aligned against the draft assembly and the spinner scaffolder looks for
links between contigs and create scaffolds. 

Original pipeline from Zemin Ning (zn1@sanger.ac.uk): http://www.sanger.ac.uk/science/tools/smis.
Version here modified to use the bwa aligner instead of smalt, and to automize
compiling and running.

### COMPILE: 
SMIS requires zlib, make sure it's in your PATH.
SMIS also requires bamtools, please create a variable MYBAMTOOLS
with the location of your bamtools:

	$ export MYBAMTOOLS=/full/path/to/bamtools

After creating MYBAMTOOLS, compile the codes by:

	$ ./makeall.sh


### RUN 
#### Step 1:	
   	
	create a variable MYSMISDIR with the location of your SMIS folder:
 	$ export MYSMISDIR=/full/path/to/smis_folder

#### Step 2:

	$MYSMISDIR/setup.sh </full/path/to/destdir> <draft_assembly> <long_reads>
	where:
   	   /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
   	   draft_assembly: fasta file of the assembly to be scaffolded
  	   long_reads: fastq file of long reads for scaffolding

#### Step 3:
   
	cd /full/path/to/destdir
   	./mysmissv.sh

### Results

Scaffolds will be in /full/path/to/destdir/spinner_scaffolds.fasta


