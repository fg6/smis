//struct _Scaffold;
//struct _Distribution;

typedef struct _PairLib
// info common to all reads in a particular library
// a pointer to such a stucture is included in each read structure
{
	int32my std;			// std dev of insert size for pair
	int32my size;			// insert size for pair
	float weight;			// weight of pair when estimating gap sizes
	int32my len;			// expected length of reads
	boolmy_fast orientation;	// 0="innie", 1="outie" 
	Distribution dist;		// distribution of insert sizes (see gapest.h).
	int32my libIndex;		// where does this lib come in the array of libs
	int32my numLibs;		// the total number of libraries used (duplicated in every lib, but anyway)
} PairLib;



typedef struct _Read
// info on read's placement on contig (and scaffold) and so on
{
	//note that the readname is not reatined here.
	char dir_c;                  //direction, +1 for forward and 0 for backwards
//	char contigName[ReadNameLen];   //name of contig from smort
	//struct _Scaffold *scaffold;     //index of scaffold, (initially scaffolds contain one contig and have the same order)
	int32my scaffoldIndex;  //index of the containing contig on the initial array of contigs.  Also of scaffold in INITIAL list.
	int32my matchStart;                 //start of match on contig
	int32my matchEnd;                   //end of match on contig
	char dir_s;			//direction relative to this programs's scaffold -- intially the same as dir_c but changes as contigs are flipped and joined
	int32my scaffoldStart;             // start of match on this programs's scaffolds
	int32my scaffoldEnd;
	char swScore;                   //Smith Waterman score
	char matchScore;                //match quality score
	char chimeric, lowQual, crossBiotin; //bool flags set for various maladies.
	struct _Read *next;		//reads are also in linked lists of reads associated to each edge.
	struct _Read *mate;		//makes things slightly faster to store the mate position rather than finding out where it is each time.
	PairLib *lib;			//data common to all reads in a library, see above.
} Read;

typedef struct //one of these will store all read data
{
        Read *read;   //to array of reads
        int64my readNum; //number of reads in it
	PairLib *lib; //array of values for data common to all reads in library.
} ReadList;

typedef struct _NameTable //a list of these gives the contig names in alpahbetical order.
{
        char name[ReadNameLen];     //contig name from fastq
        int32my index; //there will not be more than 2 billion contigs one hopes.
} NameTable;

typedef struct _Contig//stores things about contig necessary for building scaffolds
{
        char name[ReadNameLen];   //from the original smort/fastq files
        int32my len;       //number of bases in this contig
        int32my numReads;       //number of reads matching the contig -- ***this includes pairs not spanning two contigs***.
	char *bases;		//string containing base info
	char *quals;		//string containing qual info
} Contig;

typedef struct _ContigList //One of these will be the struct conataining all the contig data
{
        Contig *contig;   //position of first contig in _array_ of contigs, not linked list
        int32my contigNum;    //number of contigs in list
	char *FastQData;  // string containing fastq file (contains base and qual info)
	char FastQDataRetained; //set when the former is allocated.
	NameTable *nameTable; //a table of names along with the pointer to the appropriate scaffold.
	int32my numInName;     //if there is a unique number in the name of each contig this says where it starts in the string or -1 if not
	int32my *numToIndex; // NumToIndex[i]=j, i is the number in the contig's name and j is the index in the initial array of contigs.
	float readsPerBaseLn2; //the number of reads mapped per base, mulitplies by ln 2
} ContigList;

//external Functions
void *MallocWithAssert(int64my size, char *msg);

void *ReallocWithAssert(void *p_in, int64my size, char *msg);

Read *GetMate(Read *read); //gets pointer to mate of given read.

ReadList *LoadSmalt(char *fname, ReadList *readList, ContigList *contigList, PairLib *lib);
//reads data to readList from ONE smort file (about placement of reads on contigs and pair info).  Files should be preprocessed to contain nothing but cigar:D lines

void InitialiseReadList(ReadList *readList);
// sets number of reads to 0 and an empty list of reads

boolmy_fast IsCrossBiotin(Read *read1,Read *read2, int16f crossBiotinQualityDeficit);
// 1 if the read is judge cross biotin and 0 if not.
// a good value for crossBiotinQualityDeficit is 5. 
// the threshold readlen<=100 ? readLen - 5 : 96
// comes from Zemin's emperical values:
//"SW 70 for 2x76bp data and SW 65 for 2x70bp data. You may also set SW=95 for 2x100bp datasets."
// but for longer reads this may need to be larger I suppose.

void DeleteReadList (ReadList *readList);
//free memory