#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "gapest.h"
#include "readlist.h"
#include "scaffoldlist.h"
#include "refine.h"

// file contains the main routine and functions that process the files
// containing the parameters.

void DisplayUsage()
{
	fprintf(stderr,"Usage: spinner -f <files file name> -s <settings file name> \n");
	fprintf(stderr,"Scaffolds contigs from fastq or fasta file using smalt (bwa) output \n");
	fprintf(stderr,"Files file contains details of smalt (bwa) output, fast files and output names.\n");
	fprintf(stderr,"default names files.txt and settings.txt\n");
}

int ProcessCommandLine(int argc, char *argv[], char *filesName, char *settingsName)
{
// first a function to handle command line input
// which sets the names of the two parameter files.
	char **thisArg;
	char option_c;
	int numLeft;

	numLeft=argc;
	thisArg=argv;

	if (numLeft == 2 && ( !strcmp(*(argv+1),"--help") || !strcmp(*(argv+1),"-?") || !strcmp(*(argv+1),"-h") ) )
		{DisplayUsage(); exit(1);}

	++ thisArg; //get rid of executable name from argument list
	-- numLeft;

	strcpy(filesName,DEFAULT_FILES_NAME);
	strcpy(settingsName,DEFAULT_SETTINGS_NAME);
	while (numLeft > 0 && **thisArg == '-')
	{
		option_c=*((*thisArg) + 1);
		++ thisArg;
		-- numLeft;
		switch ( option_c ) {
		case 'f':
			strcpy(filesName,*thisArg);
			++ thisArg;
			-- numLeft;
			break;
		case 's':
			strcpy(settingsName,*thisArg);
			++ thisArg;
			-- numLeft;
			break;
		default:
			fprintf(stderr,"Unknown option: -%c .\n",option_c); DisplayUsage(); exit(1);
		}
	}
	if (numLeft > 0) {fprintf(stderr,"Too many arguments for option -%c .\n",option_c); DisplayUsage(); exit(1);}
	//printf("Files given in file: %s\n", filesName);
	//printf("Settings given in file: %s\n", settingsName);

	return (0);
}

int LoadFiles(char *fname, char *fastqOutName, char *fastaOutName, char *contigOutName, char *contigBackwardsOutName, char *GDFOutName, ReadList *readList, ContigList *contigList)
// looks at fname, file containing file names and associated details
// and loads the input files.
{
	FILE *fp;
	int i;
	char line[MaxLineLen], inputName[MaxLineLen],token[MaxLineLen], orientString[10];
	char filesName[MaxLineLen], settingsName[MaxLineLen];  //the two files containing the set-up data
	int32my insertSize, std, len; //stats for the smalt (bwa) files
	float weight;
	int32my loadedFastfile,loadedSmalt,fastaIn;
	int16f numSmalts=0;
	PairLib *currentLib, *iLib;

	fp = fopen(fname,"r");
        if (fp==NULL) {fprintf(stderr,"Could not open file %s\n",fname); exit(1);}

	loadedFastfile=0;
	loadedSmalt=0;
	fastaIn=0;
	strcpy(fastqOutName,NO_NAME);
	strcpy(fastaOutName,NO_NAME); //in this case this is ajust a mark saying that no file need be produced.
	strcpy(contigOutName,NO_NAME);
	strcpy(contigBackwardsOutName,NO_NAME);
	strcpy(GDFOutName,NO_NAME);
	contigList->numInName=6; // 7th indexed from 0

	for(i=0; fgets(line,MaxLineLen,fp) != NULL; i++){
	   //if (!strncmp(line,"bwa",3)) ++numSmalts;
           if (!strncmp(line,"bwa",3) || !strncmp(line,"smalt",5)) ++numSmalts; 
        }
	readList->lib=MallocWithAssert(sizeof(PairLib) * numSmalts,"LoadFiles, readList->lib");
	rewind(fp);



	currentLib=readList->lib;
	for(i=0; fgets(line,MaxLineLen,fp) != NULL; i++)
        {

		sscanf(line,"%s", token);
		if (!strncmp(token,"#",1)) continue;
		if (!strcmp(token,"smalt") || !strcmp(token,"bwa") ) {
                //if (!strcmp(token,"smalt (bwa)")){
			 printf("%s \n", token);
			loadedSmalt=1;
			if (!loadedFastfile) {fprintf(stderr,"in %s: smalt (bwa) file given before fastq/a \n",fname); exit(1);}
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: smalt (bwa) token but no following line\n",fname); exit(1);}
			sscanf(line,"%s %d %d %f %d %s", inputName, &insertSize, &std, &weight, &len, orientString);
			currentLib->std=std;
			currentLib->size=insertSize;
			currentLib->weight=weight;
			currentLib->len=len;
			currentLib->libIndex=currentLib - readList->lib;
			if (!strcmp(orientString,"in")) currentLib->orientation=0; else
				if (!strcmp(orientString,"out")) currentLib->orientation=1; else
					{fprintf(stderr,"in %s: smalt (bwa) data line should be <filename> <insert size> <standard deviation> <weight> <read length> <orientation, must = \"in\" | \"out\"> \n",fname); exit(1);}
			printf("Smalt (bwa) file %s with insert size %d.  Loading. \n", inputName, insertSize);
			LoadSmalt(inputName, readList, contigList, currentLib);
			printf("Total number of reads is now %d \n",readList->readNum);
			++currentLib;
		}
		if (!strcmp(token,"fastq-in")) {
			//printf("fastq-in in input\n");
			if (loadedFastfile) {fprintf(stderr,"in %s: fastq given after another fast file \n",fname); exit(1);}
			loadedFastfile=1;
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: fastq-in token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", inputName);
			printf("Loading fastq file %s\n", inputName);
			LoadFastq(inputName, contigList);
			//printf("Number of Contigs: %d \n", contigList->contigNum );
		}
		if (!strcmp(token,"fasta-in")) {
			//printf("fastq-in in input\n");
			if (loadedFastfile) {fprintf(stderr,"in %s: fasta given after another fast file \n",fname); exit(1);}
			loadedFastfile = fastaIn = 1;
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: fastq-in token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", inputName);
			printf("Loading fastq file %s :\n", inputName);
			LoadFasta(inputName, contigList);
			//printf("Number of Contigs: %d \n", contigList->contigNum );
		}
		if (!strcmp(token,"fastq-out")) {
			if (fastaIn) {fprintf(stderr,"in %s: fastq output requested but some inputs are fasta. fastq->fasta not possible. \n",fname); exit(1);}
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: fastq-out token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", fastqOutName);
		}
		if (!strcmp(token,"fasta-out")) {
			//printf("fastq-out in input\n");
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: fasta-out token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", fastaOutName);
		}
		if (!strcmp(token,"contigs-out")) {
			//printf("contigs-out in input\n");
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: contigs-out token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", contigOutName);
		}
		if (!strcmp(token,"contigs-backwards-out")) {
			//printf("contigs-out in input\n");
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: contigs-backwards-out token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", contigBackwardsOutName);
		}
		if (!strcmp(token,"GDF-out")) {
			//printf("contigs-out in input\n");
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: contigs-backwards-out token but no following line\n",fname); exit(1);}
			sscanf(line,"%s", GDFOutName);
		}
		if (!strcmp(token,"number-in-contig-name")) {
			//printf("contigs-out in input\n");
			if (fgets(line,MaxLineLen,fp)==NULL) {fprintf(stderr,"in %s: number-in-contig-name token but no following line\n",fname); exit(1);}
			sscanf(line,"%d", & (contigList->numInName));
			-- (contigList->numInName); //from normal counting to c indexing.
			printf("Will index contigs by number in name string at (c) position %d \n", contigList->numInName);
		}
	}

	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
	for (iLib=readList->lib;iLib <= currentLib; ++iLib) iLib->numLibs=currentLib - readList->lib;

	if ( contigList->numInName == 7 ) printf("Assuming standardised names for contigs <six characters><contignumber><other> \n", contigList->numInName);
	printf("Fastq file to be written to: %s\n", fastqOutName);
	printf("Fasta file to be written to: %s\n", fastaOutName);
	printf("Contig data file to be written to: %s\n", contigOutName);
}

#define BAD_PARAMETER -1
#define PAIR_THRESHOLD_CODE 0
#define STRENGTH_TEST_RATIO_CODE 1
#define EXCLUDE_SMALL_CONTIGS_CODE 2
#define PREVENT_TWO_SMALL_SCAFFOLDS_JOINING_CODE 3
#define FIND_SIGMAS_CODE 4
#define MIN_UNCERTAINTY_CODE 5
#define JUMP_LENGTH_CODE 6
#define GOOD_GAP_ESTIMATES_CODE 7
#define VERBOSE_CODE 8
#define CHIMERA_SIGMAS_CODE 9
#define MATCH_SCORE_THRESHOLD_CODE 10
#define CROSS_BIOTIN_PARAMETER_CODE 11
#define DIRECTION_TEST_RATIO_CODE 12
#define WEAK_EDGE_STRENGTH_RATIO_CODE 13
#define STAGE_CODE 14

int32my GetCode(char *token)
{
	if (!strcmp(token,"PAIR_THRESHOLD")) return(PAIR_THRESHOLD_CODE);
	if (!strcmp(token,"STRENGTH_TEST_RATIO")) return(STRENGTH_TEST_RATIO_CODE);
	if (!strcmp(token,"EXCLUDE_SMALL_CONTIGS")) return(EXCLUDE_SMALL_CONTIGS_CODE);
	if (!strcmp(token,"PREVENT_TWO_SMALL_SCAFFOLDS_JOINING")) return(PREVENT_TWO_SMALL_SCAFFOLDS_JOINING_CODE);
	if (!strcmp(token,"FIND_SIGMAS")) return(FIND_SIGMAS_CODE);
	if (!strcmp(token,"MIN_UNCERTAINTY")) return(MIN_UNCERTAINTY_CODE);
	if (!strcmp(token,"JUMP_LENGTH")) return(JUMP_LENGTH_CODE);
	if (!strcmp(token,"GOOD_GAP_ESTIMATES")) return(GOOD_GAP_ESTIMATES_CODE);
	if (!strcmp(token,"VERBOSE")) return(VERBOSE_CODE);
	if (!strcmp(token,"CHIMERA_SIGMAS")) return(CHIMERA_SIGMAS_CODE);
	if (!strcmp(token,"MATCH_SCORE_THRESHOLD")) return(MATCH_SCORE_THRESHOLD_CODE);
	if (!strcmp(token,"CROSS_BIOTIN_PARAMETER")) return(CROSS_BIOTIN_PARAMETER_CODE);
	if (!strcmp(token,"DIRECTION_TEST_RATIO")) return(DIRECTION_TEST_RATIO_CODE);
	if (!strcmp(token,"WEAK_EDGE_STRENGTH_RATIO")) return(WEAK_EDGE_STRENGTH_RATIO_CODE);
	if (!strcmp(token,"STAGE")) return(STAGE_CODE);
	return (BAD_PARAMETER);
}

void LoadSettings(char *fname, MergeParams **paramsPtr,int32my *numStages)
{
	FILE *fp;
	int32my i=0, firstStage=1;
	char token[MaxLineLen]="initialise to stop valgrind freaking out", value[MaxLineLen]="stop valgrind freaking out";
	char line[MaxLineLen]="stop valgrind freaking out";
	int32my turnOn, turnOff;
	float number, minThreshold;
	MergeParams *params;

	fp = fopen(fname,"r");
        if (fp==NULL) {fprintf(stderr,"Could not open file %s\n",fname); exit(1);}

	*numStages=0;
	while(fgets(line,MaxLineLen,fp) != NULL)
		if (!strncmp(line,"STAGE",5)) ++ *numStages;
	rewind(fp);

	*paramsPtr=(MergeParams *) MallocWithAssert(sizeof(MergeParams) * *numStages,"LoadSettings, params");
	params=*paramsPtr;

	while(fgets(line,MaxLineLen,fp) != NULL)
        {
		if (*line=='\n' || *line=='#') continue;
		sscanf(line,"%s %s", token, value);
		if (!strcmp(value,"ON"))  turnOn=1; else turnOn=0;
		if (!strcmp(value,"OFF")) turnOff=1; else turnOff=0;
		number=atof(value);

		switch (GetCode(token)) {
		case STAGE_CODE:
			if (firstStage) {
				//fill in default values
				params[0].edgeParams.matchScoreThreshold=15;
				params[0].edgeParams.chimeraSigmas=5;
				params[0].edgeParams.numPairThreshold=30;
				params[0].edgeParams.StrengthRatioForWeakEdges=-1.0;
				params[0].edgeParams.crossBiotinSWDeficit=20;
				params[0].findSigmas=5.0;
				params[0].strengthRatio=0.0 ;
				params[0].minUncertainty=100;
				params[0].excludeSmallContigs=0;
				params[0].smallContigThreshold=3000;
				params[0].preventTwoSmallScaffoldMergeThreshold=0;
				params[0].dirRatio=3.0;
				params[0].skipLen=NO_SKIP;
				params[0].mergeMark=1;
				params[0].maxScaffoldingSteps=20;
				params[0].minMergesPerStep=1;
				params[0].verbose=0;
				params[0].useMLE=1;
				firstStage=0;
			} else {
				++i;
				params[i]=params[i-1]; //copy previous settings forward
				params[i].mergeMark=i+1;
			}
			break;
		case PAIR_THRESHOLD_CODE:
			params[i].edgeParams.numPairThreshold=(int32my) number;
			break;
		case STRENGTH_TEST_RATIO_CODE:
			params[i].strengthRatio=number;
			break;
		case EXCLUDE_SMALL_CONTIGS_CODE:
			if (turnOff) params[i].excludeSmallContigs=0; else {
				params[i].excludeSmallContigs=1;
				params[i].smallContigThreshold=(int32my) number;
			}
			break;
		case PREVENT_TWO_SMALL_SCAFFOLDS_JOINING_CODE:
			if (turnOff) params[i].preventTwoSmallScaffoldMergeThreshold=0;
			else params[i].preventTwoSmallScaffoldMergeThreshold=(int32my) number;
			break;
		case FIND_SIGMAS_CODE:
			params[i].findSigmas=number;
			break;
		case MIN_UNCERTAINTY_CODE:
			params[i].minUncertainty=(int32my) number;
			break;
		case JUMP_LENGTH_CODE:
			if (turnOff) params[i].skipLen=NO_SKIP;
			else params[i].skipLen=(int32my) number;
			break;
		case GOOD_GAP_ESTIMATES_CODE:
			if (turnOff) params[i].useMLE=0;
			if (turnOn) params[i].useMLE=1;
			break;
		case VERBOSE_CODE:
			if (turnOff) params[i].verbose=0;
			if (turnOn) params[i].verbose=1;
			break;
		case CHIMERA_SIGMAS_CODE:
			params[i].edgeParams.chimeraSigmas=number;
			break;
		case MATCH_SCORE_THRESHOLD_CODE:
			params[i].edgeParams.matchScoreThreshold=(int32my) number;
			break;
		case CROSS_BIOTIN_PARAMETER_CODE:
			params[i].edgeParams.crossBiotinSWDeficit=(int32my) number;
			break;
		case DIRECTION_TEST_RATIO_CODE:
			params[i].dirRatio= number;
			break;
		case WEAK_EDGE_STRENGTH_RATIO_CODE:
			params[i].edgeParams.StrengthRatioForWeakEdges= number;
			break;
		default:
			printf("In settings file %s unknown parameter %s\n",fname,token);
			fprintf(stderr,"In settings file %s unknown parameter %s\n",fname,token);
			exit(1);
		}
	}
	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
	//StrengthRatioForWeakEdges must be calculated.
	minThreshold=10000;
	for(i=0;i<*numStages;++i)
		if(params[i].edgeParams.numPairThreshold < minThreshold) minThreshold= params[i].edgeParams.numPairThreshold;
	for(i=0;i<*numStages;++i)
		if(params[i].edgeParams.StrengthRatioForWeakEdges < 0.0)
			params[i].edgeParams.StrengthRatioForWeakEdges= 
				params[i].strengthRatio < minThreshold / params[i].edgeParams.numPairThreshold ?
				minThreshold / params[i].edgeParams.numPairThreshold : params[i].strengthRatio;

}

//////////////////////////////////////// main /////////////////////////////////////////

int main(int argc, char *argv[])
{
	ReadList readList;
	ContigList contigList;
	ScaffoldList scaffoldList;

	char fastqOutName[MaxLineLen], fastaOutName[MaxLineLen], contigOutName[MaxLineLen], contigBackwardsOutName[MaxLineLen]; //files to write to
	char GDFOutName[MaxLineLen];
	char filesName[MaxLineLen], settingsName[MaxLineLen]; //files containing the run configuration
	int32my i;
	int64my totalLen;
	MergeParams *params;
	int32my numStages;

	printf("++++++++++ Spinner 1.0 ++++++++++\n\n");

	// Set the structures to contain 0 reads and 0 contigs.
	InitialiseReadList(&readList);
	InitialiseContigList(&contigList);

	// Now we need to get the parameters for the run
	// And fill in the reads and contigs using the files provided.
	ProcessCommandLine(argc, argv, filesName, settingsName); 
	printf("Loading settings from file %s .\n",settingsName);
	LoadSettings(settingsName,&params,&numStages);
	printf("Loading file details from file %s .\n",filesName);
	LoadFiles(filesName, fastqOutName, fastaOutName, contigOutName, contigBackwardsOutName,GDFOutName, &readList, &contigList); //now the reads and contigs should be loaded.

	//the we create the scaffolds to each contain one contig
	printf("Creating scaffold list\n");
	SetScaffoldListFromContigList(&contigList,&scaffoldList);
	printf("Creating edges\n");
	//and firther fill in the ScaffoldList structure to define the edges.
	MakeEdgesFromReads(&readList, &contigList, &scaffoldList);
	ReportOnEdges(&scaffoldList);

	totalLen=0;
	for (i=0;i<scaffoldList.scaffoldNum;++i) totalLen+=scaffoldList.scaffold[i].len;
	printf("Aggregate contig len %" PRI64 " \n",totalLen);

	// the main procedure
	DoScaffolding(&scaffoldList, params, numStages);
	SortScaffoldListByLength(&scaffoldList);
	printf("\n");
	////////////////////////
	//printf("(VALIDATE) checking halfedge list \n \n \n \n");
	//ValidateAllHalfEdges(&scaffoldList);
	///////////////////////

	//save and deallocate
	printf("Saving: \n");
	if (strcmp(fastqOutName,NO_NAME)) SaveFastq(&scaffoldList, fastqOutName);
	if (strcmp(fastaOutName,NO_NAME)) SaveFasta(&scaffoldList, fastaOutName); 
	if (strcmp(contigOutName,NO_NAME)) SaveContigsInScaffoldsProperFormat(&scaffoldList, contigOutName, 1);	
	if (strcmp(contigBackwardsOutName,NO_NAME)) SaveContigsInScaffoldsProperFormat(&scaffoldList, contigBackwardsOutName, 0);
	if (strcmp(GDFOutName,NO_NAME)) SaveGDF(&scaffoldList, GDFOutName, 500, -1, -1);


	printf("freeing memory: \n");
	DeleteReadList(&readList);
	DeleteContigList(&contigList);
	DeleteScaffoldList(&scaffoldList);
	free(params);
	printf("done freeing memory.  End. \n");
  return EXIT_SUCCESS;
}
