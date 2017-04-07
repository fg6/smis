//lots of structs

struct _Scaffold;

typedef struct _SetEdgeParams //a list of these gives the contig names in alpahbetical order.
{
        char matchScoreThreshold; //reads with a match quality score below this are not used.
	int32my chimeraSigmas; //reads that imply a gap size more nagtive than -  ChimaraSigmas * the read insert size std dev are ignored at present.
	int16my numPairThreshold; //edges with less than this number of pairs suopporting them are marked bad quality.
	int16f crossBiotinSWDeficit; //the SW score of one read has to be this much less than the read length to trip the CrossBiotin check.
	float StrengthRatioForWeakEdges; //similar to StrengthRatio for edges that fall below the strength threshold but can still block other edges.
} SetEdgeParams;

struct _HalfEdge;

typedef struct _ContigInScaffold //doubly-linked lists representing the contigs making up a scaffold 
{
        Contig *contig;   // contig that's at this position in Scaffold
	boolmy_fast dir;          // direction of contig compared to its orginal orientation (0=reversed, 1=normal)
        struct _ContigInScaffold *forwards;    //next contig TOWARDS front (orginally to the RIGHT) or null if front
        struct _ContigInScaffold *backwards;    //next towards back or null if back
	int32my forwardGap; //estimated size of gap to the next contig forward or 0 if front
	int32my backwardGap; //estimated size of gap to the next contig backward or 0 if back
	char backwardMark, forwardMark; //used to mark the successful merges with some label, now the merge stage at which they ocurred, which appears in some output files so that you can check what happenned and when.
	struct _HalfEdge *forwardHe, *backwardHe; //the halfEdges that have been merged on or NULL at ends
} ContigInScaffold;

typedef struct _Edge// labels for the full edge (information that is shared by the two half edges)
{
	//int isNotAltered; //true if the set has been set, and has not been made in need of setting again (i.e. merged with some other edge) since.
	int16f numPair;			// number of pairs confirming the link, including all types of bad link
	int16f numGoodPair;		// number of pairs confirming the link including no bad links
	int16f numBiotin; // number of cross biotin read pairs
	boolmy_fast edgeQual; // after some tests 0=suspect 1=okay, at present just a threshold on number of good pairs
	boolmy_fast conflictCheck; //if this edge should be condiered when calculating conflicts, just a lower threshold.
	//NB: the idea here is that, if the threshold is 30, an edge of strength 29 shouldn't be merged but could still block an edge of strength 30
	boolmy_fast isMultiEdge; //set if there is another edge between the same two scaffolds with different directions on it.
	float gapSize;			// esimate of gap size
	float w;				// inverse of std dev of this gap size iteratively built up. 
	float meanScore;
} Edge;

typedef struct _HalfEdge  //a half-edge of the bidirectional edge, to be in a list in the scaffolf associated to it
{
	struct _Scaffold *oppositeScaffold;	//at the other end of the edge
	struct _HalfEdge *otherEnd;		//the corresponding halfedge from that scaffold's vertex
	boolmy_fast dir;				//whether half-edge enters front or back of this scaffold
	int numReads;			//number of supporting reads
	Edge *edge;				//points to the associated edge
	Read *read;				//list of reads associated to the half edge (the relevant halves of the mate pairs).
	struct _HalfEdge *next;			//next in scaffold's list of halfedges
	struct _HalfEdge *prev;			//previous in scaffold's list of halfedges
} HalfEdge;

typedef struct _Scaffold // linked list elements for the scaffolds, or vertices of graph
{
	//The data for this scaffold as vertex of the graph
	int frontEdgeNum;		//number of edges to this scaffold's front end
	int backEdgeNum;		//number of edges to this scaffold's back end
	HalfEdge* frontHalfEdge;	//first in the linked list of half-edges into front
	HalfEdge* backHalfEdge;		//first in the linked list of half-edges into back
	//properties of scaffold
	int32my len;                          //estimated number of bases in this scaffold
	int16my contigNum;                    //number of contigs in this scaffold
	int32my numReads;       //number of reads matching scaffold -- ***this includes pairs not spanning two contigs***.
	ContigInScaffold *frontContig;   //points to struct that points to contig at "front" of list of contigs in this scaffold(orginally the rightmost).
	ContigInScaffold *backContig;    //same but points to struct that points to contig at "back" 
	struct _Scaffold *next;           //points to next scaffold in linked list of scaffolds (I require no special order, though can be sorted largest first)
	struct _Scaffold *prev;           //profiling suggested that I should make this a doubly linked list, surprising how muh time it saves.
	//data for merging process
	HalfEdge *frontClosest,*backClosest; //when the nearest scaffolds to this one are calculated the edges to them are recorded here.
	boolmy isDormant ; // turned to 1 (default 0) if we don't want to merge his scaffold (if it is a repeat we want to jump over for instance)
	char mark;			//used now to mark connected components for saving the graph as GDF
	boolmy isFrontJunction;	// true if there is an unresolvable conflict between connections to the front
	boolmy isBackJunction;	// similar
} Scaffold;

typedef struct //all this gives this list of scaffolds (which is the graph basically).
{
	Scaffold *scaffold;             // position of first contig in linked list, intially the same as...
	Scaffold *firstScaffold;             // position of first contig in memory block
	int32my scaffoldNum;                        // number of scaffolds in the list INITIALLY
	int64my edgeNum;                        // number of edges in the list INITIALLY
	//the scaffold list elements contain lists of contigs and edges which must be memory managed.
	ContigInScaffold  *firstContigInScaffold;				// block of memory where all the ContigInScaffold structs are actually located.
	ContigInScaffold  *newContigInScaffold; // next empty location for one (more of these are created as we go along).
	HalfEdge *firstHalfEdge; //block for half edges (so we can deallocate it later)
	Edge *firstEdge; //block for edges.
} ScaffoldList;


//*********External Functions**********************************//

void CrossBiotinFlag(HalfEdge *he, int16f SWdeficit);
//checks if reads supporting the half edge are corssBiotin and flags them accordingly

void LowMatchScoreFlag(HalfEdge *he, char threshold);
//checks if reads supporting the half edge have low match qualtity score and flags them accordingly

//int HeListGood(HalfEdge *he)

HalfEdge *GetHalfEdge(Scaffold *scaf1, boolmy_fast dir1, Scaffold *scaf2, boolmy_fast dir2);
//  if there is a half edge from scaf1 to scaf2 with given directions, returns a pointer to it, and NULL if not.

void SortScaffoldListByLength(ScaffoldList *scaffoldList);
// sorts the linked list of scaffolds so that the largest is the head and so on.

void MakeSmallContigsDormant(ScaffoldList *scaffoldList, int smallContigThreshold);
// turns on the "dormant" flag for all small contigs.

void MakeJunctionsDormant(ScaffoldList *scaffoldList);

void ResetDormant(ScaffoldList *scaffoldList);

void SetEdge(HalfEdge *halfEdge,Edge *edge, SetEdgeParams *params);
//sets the gapSize and all other data in edge, based on read list in halfedge.

void QuickEstimateGapNoCutoff(HalfEdge *he);
//finds the gap for the relevant half edgeby weighted average without throwing away "chimeric" edges

void RemoveEdge(HalfEdge *he1);
//remove the two halfedges (one of which is he1)
// from the lists of edges associated to the two scaffolds
//has to be in properly constructed scaffold list or something odd will happen

void RemoveEdgeSafe(HalfEdge *he1, Scaffold *scaf1, Scaffold *scaf2);
// remove the two halfedges (one of which is he1)
// from the lists of edges associated to the two scaffolds
// works when the oppositeScaffolds have been made incorrect
// But of course the scaffolds must be correct.

int RemoveScaffold(ScaffoldList* scaffoldList, Scaffold *scaf);
//remove the scaffold from the scaffoldlist
//has to be in the list, or it returns 0

Contig *GetContigFromName(char *name, ContigList *contigList);
//returns pointer to contig of the given name.

Scaffold *GetScaffold(Read *read1, ContigList *contigList, ScaffoldList *scaffoldList);
//supplies pointer to the scaffold on which read 1 is matched
//returns NULL if the name is not in the list (this happens sometimes as some contigs are not good and are not included in fastq file).

void SetScaffoldIndex(Read *read, char *contigName, ContigList *contigList);

ContigList *LoadContigListFromFastq(char *fname, ContigList* contiglist);
//uses fastq file purely to get lengths of the contigs,
//the only thing we need at this stage that is not in smort file

boolmy_fast LoadFastq(char *fname, ContigList *contigList);
// fills in config info from fastQ file: names, lengths and base and qual info
// allocates memory block pointed to by contigList->fastQData

boolmy_fast LoadFasta(char *fname, ContigList *contigList);
// fills in config info from fasta file: names, lengths and base info
// allocates memory block pointed to by contigList->fastQData

void InitialiseContigList(ContigList *contigList);
// sets num of entries in list to zero and makes list empty

ScaffoldList *SetScaffoldListFromContigList(ContigList* contiglist,ScaffoldList* scaffoldlist);
//intialises contigList, creating one scaffold for each contig containing only that contig

int SaveFastq(ScaffoldList *scaffoldList, char *fname);

int SaveFasta(ScaffoldList *scaffoldList, char *fname);

int SaveGDF(ScaffoldList *scaffoldList, char *fname, int min_size, int num, int depth);
// saves a graph file showing the scaffolding graph
//if num==0 saves all, if > 0 saves num connected components only
//scaffolds of length < min_size are ignored.

int SaveContigsInScaffoldsProperFormat(ScaffoldList *scaffoldList, char *fname, int dir);

ScaffoldList *MakeEdgesFromReads(ReadList *readList, ContigList *contigList,ScaffoldList *scaffoldList);
// The ScaffoldList should already have been intialised (with each scaffold containing one contig and no edges).
// Now the edges, relations between the scaffolds inferred from the read pairs, are added to the structure to make the bidirectional graph.
// Using a soft screen approach, all read pairs are recorded in the graph, even onces with bad scores.
// However, the number of pairs supporting an edge with good scores is recorded in the edge.

void DisplayClosestScaffolds(ScaffoldList *scaffoldList);
//prints out a list of scaffolds together with their front and back closest and directional properties

int DeleteScaffoldList(ScaffoldList *scaffoldList);
//free memory in scaffoldList

int DeleteContigList(ContigList *contigList);
//free memory in contigList

float EstimateRelativePosition(Scaffold *scaf, HalfEdge *he);

boolmy IsRoomInGap(Scaffold *scaf, ContigInScaffold **cis, int32my *posFrontOfGap, int32my range_min, int32my range_max, int32my size);