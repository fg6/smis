// readlist.c
// handles things specific to the read and readlist structures
// such as loading from smort files.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "gapest.h"
#include "readlist.h"
#include "scaffoldlist.h"


void *MallocWithAssert(int64my size, char *msg)
{
	void *p;
	p=malloc(size);
	if (p==NULL) {fprintf(stderr, "malloc failed: %s\n", msg); exit(1);}
	return (p);
}

void *ReallocWithAssert(void *p_in, int64my size, char *msg)
{
	void *p;
	p=realloc(p_in, size);
	if (p==NULL) {fprintf(stderr, "realloc failed: %s\n", msg); exit(1);}
	return (p);
}

Read *GetMate(Read *read)
{
	return (read->mate);
	//return ( (read - readList->read) % 2 == 0 ? read+1 : read-1); //read pairs are stored together as pairs in memory
}

boolmy_fast IsCrossBiotin(Read *read1,Read *read2, int16f crossBiotinSWDeficit)
{
// 1 if the pair is judged cross biotin and 0 if not.
// judged bad if one read has good SW score and the other is below readLen - crossBiotinSWDeficit
// 5 seems to be a good value.
//"SW 70 for 2x76bp data and SW 65 for 2x70bp data. You may also set SW=95 for 2x100bp datasets."
//may need to be revised for longer reads.
	int32f len;
	if(crossBiotinSWDeficit==-1 || (len=read1->lib->len) > 500) return 0; // if its above 500 this isn't data for which I know there are cross biotins like this 
	if (read1->swScore >= len-3 && read2->swScore < len - crossBiotinSWDeficit ) return 1;
	if (read2->swScore >= len-3 && read1->swScore < len - crossBiotinSWDeficit ) return 1; else return 0;
}

void InitialiseReadList(ReadList *readList)
{
// sets number of reads to 0 and an empty list of reads
	readList->readNum=0;
	readList->read=NULL;
}

ReadList *LoadSmalt(char *fname, ReadList *readList, ContigList *contigList, PairLib *lib)
//reads data to readList from a smalt file (about placement of reads on contigs and pair info).  Files should be preprocessed to contain nothing but cigar:D lines.
//NB: nothing in this code relies on the form of read or contig names.  But instead it does rely on readpairs being together in the file.
{
        FILE *fp;
        char line[MaxLineLen],cigar[12],readDir[2],contigName[ReadNameLen]; //halfway houses for file data.
	
        int64my i;
        Read *read;
	Contig *contig;
	//long readStart, readEnd; //because I'm not sure how to ensure that scanf
	int64my oldNum, readSum, lenSum, currentLine;
	int32my index;
	int32my firstOfPairRead;
	float toDisplay=0.0;
        //---open file------------------------

        fp = fopen(fname,"r");
        if (fp==NULL) {fprintf(stderr,"Could not open file %s\n",fname); exit(1);}

        //---find # of lines=# of reads in inter-contig pairs in file and add to # reads------------------------
        i = 0;
        while(fgets(line,MaxLineLen,fp)!=NULL) if (!strncmp(line,"cigar:D",7)) ++i;
        rewind(fp);

	oldNum=readList->readNum;
        readList->readNum += i; //stores it where it can be got at later
        //---allocate mem for readList------------------------
        if (readList->read!=NULL)
		readList->read  = (Read *) ReallocWithAssert(readList->read, sizeof(Read) * readList->readNum,"LoadSmalt, readList->read");
	else
		readList->read  = (Read *) MallocWithAssert(sizeof(Read) * readList->readNum,"LoadSmalt, readList->read");
	printf("----- LoadSmalt: allocated %" PRI64 " bytes for %" PRI64 " reads\n",(int64my) sizeof(Read) * readList->readNum,readList->readNum);

        //---go back through file and fill in ReadList------------------------
	i=oldNum;
	firstOfPairRead=0;
        for(read=readList->read + oldNum, currentLine=1; fgets(line,MaxLineLen,fp) != NULL; ++currentLine ) {
		if (!strncmp(line,"cigar",5)) { //other lines are just junk e.g. smalt output stats
			//cigar:D pairs join two contigs, A B C and S are mapped reads that don't do this,
			// N is not mapped (ignored here).
			if (line[6]=='A' || line[6]=='B' || line[6]=='C' || line[6]=='S') {
				//sscanf(line,"%*s %*s %*d %*d %*s %s %*" SCN32 " %*" SCN32 " %*s %*d", contigName );
				sscanf(line,"%*s %*s %*s %*s %*s %s %*s %*s %*s %*d", contigName );
				index=GetIndexFromName(contigName, contigList);
				if (index!=-1) ++ contigList->contig[index].numReads; //this number counts ALL mapped reads for the A-stat
			} else if (line[6]=='D') {
				firstOfPairRead = 1 - firstOfPairRead;
		                //sscanf(line,"%s %*s %*d %*d %s %s %" SCN32 " %" SCN32 " %*s %d",
				sscanf(line,"%s %*s %*s %*s %s %s %" SCN32 " %" SCN32 " %*s %d",
					cigar,
					//read[i].readName,
					//&readStart,
					//&readEnd, //no use for the start and end of read match data now
					readDir,
					contigName,
					&read->matchStart,
					&read->matchEnd,
					&read->swScore
					);
                		read->matchScore=atoi(cigar+8); //first token eg "cigar:D:59" 59 is score
				read->lib=lib;
				if (lib->orientation) { //1 here is "outie" (nonstandard)
					read->dir_c=read->dir_s= *readDir == '+' ? 1 : 0 ;
					read->scaffoldStart=read->matchEnd;
                			read->scaffoldEnd=read->matchStart;
				} else { // "innie", normal case, from which I take the terminology
					read->dir_c=read->dir_s= *readDir == '+' ? 0 : 1 ;
					read->scaffoldStart=read->matchStart;
                			read->scaffoldEnd=read->matchEnd;
				}
				read->chimeric=read->lowQual=read->crossBiotin=0;
				index=GetIndexFromName(contigName, contigList);
				//printf("contig name %s index %" PRI32 "\n", contigName, index);
				//if (! firstOfPairRead ) printf("\n");
				if (index >= contigList->contigNum) index=-1;
				if (index!=-1) ++ contigList->contig[index].numReads;
				read->scaffoldIndex=index;
				/*if ( (firstOfPairRead && i % 2) ||  (! firstOfPairRead && !(i % 2)) ) {
					fprintf(stderr,"Confused pair count at line %" PRI64 " of %s\n",currentLine,fname);
					exit(1);}*/
				read->mate= firstOfPairRead  ?  read + 1 : read - 1 ;
				if ( (! firstOfPairRead) && read->scaffoldIndex!=-1 && read->scaffoldIndex ==  read->mate->scaffoldIndex) {
					fprintf(stderr,"Cigar D pair with  same scaffold index at line  %" PRI64 " of %s\n",currentLine,fname);
					fprintf(stderr,"Probable cause: nonstandard names (not 6 letters and then unique number)\nEither rename-contigs or set number-in-contig-name in files.txt to \ncharacter at which unique number appears in contig name or -1\n");
					fprintf(stderr,"indices  %" PRI32 " and %" PRI32 "\n",read->mate->scaffoldIndex,read->scaffoldIndex);
					fprintf(stderr,line);
					exit(1);
				}
				++read; //iterate pointer and counter for useful D reads
				++i;
        		}
		}
		if (line[6]!='D' && firstOfPairRead) {
			fprintf(stderr,"Paired cigar:D entry not present line %" PRI64 " of %s\n",currentLine,fname);
			exit(1);}
	}

	if (i==oldNum) {	
		fprintf(stderr,"Smalt output file %s has added no useful read pairs \n",fname);
		fprintf(stderr,"Probable causes: names of contigs in fasta/fastq file and smalt output\n");
		fprintf(stderr,"do not match(renamed?); no cigar:D lines in file corresponding \n");
		fprintf(stderr,"to given contigs for some other reason.\n");
		fprintf(stderr,"Check names match, or try removing this file.\n");
		exit(1);
	}

	for(read=readList->read, i=0; i<readList->readNum; ++read, ++i ) //in case the memory block moved.
		read->mate= i%2 ? read-1 : read+1;

	if (firstOfPairRead) { fprintf(stderr,"half of D pair not present ln %" PRI64 " %s\n",currentLine,fname); exit(1);}
	printf("----- LoadSmalt: read %" PRI64 " lines, %" PRI64 " cigar:D lines \n", currentLine, i- oldNum);
	//now the mate is a property of the read, saving a large number of comparisions later.
	//for(read=readList->read, i=0; i < readList->readNum; ++read, ++i ) read->mate= i % 2  ?  read - 1 : read + 1 ;
	//ln2=0.69314718;
	lenSum=0; readSum=0;
	for(contig=contigList->contig, i=0; i < contigList->contigNum; ++contig, ++i ) {lenSum+=contig->len; readSum+=contig->numReads;}
	contigList->readsPerBaseLn2= toDisplay= (readSum / (float) lenSum) / 0.69314718; //ln 2.
	printf("----- LoadSmalt: reads per base / ln2 = %f\n", toDisplay);

        //---close file------------------------
        
        if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};

	//inefficient way to get the distribution, should really incorporate this into the above
	lib->dist=GetEmptyDistribution(lib->size * 2);
	if ( LoadInsertSizeSamples(&(lib->dist), fname, contigList, lib->orientation ? 'C' : 'A') < 10000) {
		//printf("Insufficient read pairs to form insert size distribution (low coverage or contigs small in comparison to insert)\n");
		printf("Substituting normal distribution for empirical distibution\n");
		SetNormalDistribution(&(lib->dist),lib->size,lib->std);
	}
	//for (i=2000; i<lib->dist.size && i<11000; ++i)
	//	printf("s %f ",lib->dist.weight[i]);
	//printf("\n");
        return(readList);
}

void DeleteReadList (ReadList *readList)
//free memory
{
	free(readList->read);
	free(readList->lib);
}

/*************************************************************************************************************************/
////////////////////////////
