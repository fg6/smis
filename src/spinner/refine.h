typedef struct _MergeParams//one of these will store all read data
{
	boolmy verbose;
	float findSigmas;  //the multiple to use as fudge factor when deciding if there is a conflict in relative positions
	float strengthRatio; //SSPACE style conflict resolution: if the ratio of the strength of two edges is smaller, resolve
	int32my minUncertainty; // minimum uncertainty to use when looking for conflict in relative positions
	boolmy excludeSmallContigs; //boolean
	int32my smallContigThreshold; //smaller than this are made dormant in this stage if excludeSmallContigs
	int32my preventTwoSmallScaffoldMergeThreshold; //if two scaffolds both of this size or smaller are found to merge, the merge is skipped.
	//(scaffold pairs of this type often have bad gap estimates and are more likely to have bad stats)
	float dirRatio;
	int32my skipLen; // when no NO_SKIP, specifies a length around contigs in which contigs and conflicts are ignored. See FindClosest.
	char mergeMark; //merges are marked with this value in the contig in scaffold list. Later you can easily see when the merge happenned from this.
	SetEdgeParams edgeParams;
	int32my maxScaffoldingSteps; //max number of times CollapseEdges should be run
	int32my minMergesPerStep; //if there are less merges than this on a step, the present scaffolding stage ends.
	boolmy useMLE; //do nice estimations of gap sizes 
} MergeParams;

//////////////////////////////////////////////////////////////////////////////
//external functions

void ValidateAllHalfEdges(ScaffoldList *scaffoldList);
//debugging: find errors in halfedge lists

void ResetEdges(ScaffoldList* scaffoldList,SetEdgeParams *params);
// each edge is marked if it is bad quality edgeQual==0
// this mark is changed to 0 for various reasons in the stage of finding closest scaffolds
// this functions sets them all back to depend only on the orginal criterion in SetEdge.

void CollapseEdges(ScaffoldList *scaffoldList, MergeParams *params);
// builds scaffolds by collapsing unambiguous edges, leaving only edges to possible junctions
// if excludeJunctions is true, scaffolds that were previously marked as junctions are not
// included when scaffolds are considered for merging

void DoScaffolding(ScaffoldList *scaffoldList, MergeParams *params, int numStages);
// the scaffolding process can be repeated with different parameters.
// this routine calls the scaffolding routine multiple times with different parameter sets.

void ReportOnEdges(ScaffoldList* scaffoldList);
////////////////////////////////////////////////////////////////////////////////////

int InsertScaffold(ScaffoldList* scaffoldList, Scaffold *bigScaf, Scaffold *littleScaf, HalfEdge *he1, HalfEdge *he2, ContigInScaffold *backwardInsertCis, int32my gapPos, MergeParams *params);

int32my GapMLE(Scaffold *scaf, HalfEdge *he);