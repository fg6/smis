
#### TO COMPILE #####
SMIS requires zlib, make sure it's in your path.
SMIS also requires bamtools, please create a variable MYBAMTOOLS
with the location of your bamtools:
$ export MYBAMTOOLS=/full/path/to/bamtools

After creating MYBAMTOOLS, compile the codes by:
   $ ./makeall.sh

The executables will be in the smissv-bin folder


#### TO RUN ####
Step 1:	
	create a variable MYSMISDIR with the location of your SMIS folder:
 	$ export MYSMISDIR=/full/path/to/smis_folder

Step 2:
	$MYSMISDIR/setup.sh </full/path/to/destdir> <draft_assembly> <long_reads>

	where:
   	   /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
   	   draft_assembly: fasta file of the assembly to be scaffolded
  	   long_reads: fastq file of long reads for scaffolding

Step 3:
   	cd /full/path/to/destdir
   	./mysmissv.sh

Scaffolds will be in /full/path/to/destdir/spinner_scaffolds.fasta


