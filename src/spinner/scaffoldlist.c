//Joe Henson
//This file deals with the Contig, ContigInScaffold, ContigList, Scaffold and ScaffoldList structures
//including Loading Contigs from fastq file and organising the data into the bidiectional graph in a ScaffoldList.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "gapest.h"
#include "readlist.h"
#include "scaffoldlist.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
/* Internal functions */

int CompareScaffoldsByLength(const void *scaf1, const void *scaf2)
//sorts greatest length first
{
	return( -( ( * (Scaffold**) scaf1)->len - ( * (Scaffold**) scaf2)->len ) ); 
}

int CompareNameTableEntries(const void *x1, const void *x2)
//for binary searches and sort of the name table, compares names alphabetically
{
	return(strcmp( ((NameTable *) x1)->name, ((NameTable *) x2)->name ) );
}

void AddEdge(HalfEdge *halfEdge, Edge **newEdgePtr)
{
// used in MakeEdgesFromReads.
// There's a block of memory allocated for Edge structures; every halfEdge created should have an Edge associated to it.
// This function takes the next empty slot from the block, and adds pointers to it from the given halfEdge and its otherEnd.
// iterates the next empty slot. (SetEdge sets the data in the edge later on).
	Edge *newEdge;
	newEdge=*newEdgePtr; // newEdge is a pointer as well actually
	++(*newEdgePtr); //this iterates the newEdge position back in the orginal pointer location.
	halfEdge->otherEnd->edge = halfEdge->edge = newEdge;	//add pointer to the new edge to the 2 half edges.	
}

void RemoveDuplicatePairs(ScaffoldList* scaffoldList)
{
// removes duplicate mate pairs from data in graph, once graph is created.
// goes through list of half edges, which contain lists of read pairs with matching scaffolds.
// So duplicates must always be found together in one half edge's list.
// takes duplicates out of these lists.

	HalfEdge *firstHalfEdge,*halfEdge;
	Read *read1, *read2, *prevRead1, *prevRead2;
	int64my numHalfEdge;
	int64my i,j,k,count;
	int32my *mStartList1,*mStartList2; //lists of the start positions of the read pairs and the number in list.
	int32my limitOnPairs=200000;
	int32my totalDuplicates;
	int32my mStart1,mStart2;

	firstHalfEdge=scaffoldList->firstHalfEdge;
	numHalfEdge=scaffoldList->edgeNum * 2;

	mStartList1= (int32my *) MallocWithAssert(sizeof(int32my) * limitOnPairs,"RemoveDuplicatePairs, mStartList1");
	mStartList2= (int32my *) MallocWithAssert(sizeof(int32my) * limitOnPairs,"RemoveDuplicatePairs, mStartList2");

	totalDuplicates=0;
	for (i=0;i<numHalfEdge; i+=2)  //halfedges are next to their other ends in the list at this stage.  Only need one of pair.
	{ //go through the list of all halfedges
		halfEdge=firstHalfEdge+i;
		count=0;
		for (read1=halfEdge->read;read1!=NULL;read1=read1->next) //list of reads in matching half edges in same order
		{ //go through the reads (and get read pairs) in the half edge's list.
			read2=GetMate(read1);
			mStart1=read1->matchStart;
			mStart2=read2->matchStart;
			for (j=0;j<count;j++) //check pair against list of non-duplicate pairs we have compiled so far
			{
				if ( (mStartList1[j] - 3 <= mStart1  &&  mStart1 <= mStartList1[j] + 3) && 
					 (mStartList2[j] - 3 <= mStart2  &&  mStart2 <= mStartList2[j] + 3)) break; //break if the approx start positions of the two reads are already in list
			}
			if (j==count)  //if we got to the end of the loop this is not a duplicate of one in the array already
			{ //so add to list and move prevRead forward for next loop
				mStartList1[count]=mStart1;
				mStartList2[count]=mStart2;
				prevRead1=read1; prevRead2=read2;
				++count;
				if (count>limitOnPairs) {
					limitOnPairs*=2;
					mStartList1= (int32my *) ReallocWithAssert(mStartList1, sizeof(int32my) * limitOnPairs,"RemoveDuplicatePairs, mStartList1");
					mStartList2= (int32my *) ReallocWithAssert(mStartList2, sizeof(int32my) * limitOnPairs,"RemoveDuplicatePairs, mStartList2");
				}
				//printf("count: %d \n",count);
			} else { //otherwise this is duplicate and we remove it from the half edge's list by joining the previous entry with the next one
				//printf("Removing duplicate pair \n"),
				prevRead1->next=read1->next; //prevRead is set at end of first exectution of loop and this condition is never satisfied on that loop
				prevRead2->next=read2->next;
				--(halfEdge->numReads);
				--(halfEdge->otherEnd->numReads);
				++ totalDuplicates;
			}
		}
		//for (read1=halfEdge->read, k=0; read1!=NULL; read1=read1->next, ++k) {}
		//if (k != halfEdge->numReads) printf("DUPLICATE FAULT %p Numreads recorded %d in list %d \n",halfEdge,k,halfEdge->numReads);
	}
	printf("(DUPLICATES) %ld duplicate read pairs found \n", totalDuplicates);
}

char ComplementBase(char base)
{
	switch( base ){
		case 'A':
			return('T');
		case 'C':
			return('G');
		case 'G':
			return('C');
		case 'T':
			return('A');
		case 'a':
			return('t');
		case 'c':
			return('g');
		case 'g':
			return('c');
		case 't':
			return('a');
		case 'N':
			return('N');
		case 'n':
			return('n');
		default:
			printf("Letter %c not known. Unchanged.\n", base);
	}
}

boolmy IsRoomInGap(Scaffold *scaf, ContigInScaffold **cis, int32my *posFrontOfGap, int32my range_min, int32my range_max, int32my size)
//if one and only one gap in the range of positions (relative to front of scaf) has size >= size returns true.
//in this case also returns a pointer to the contig in scaffold in *cis:
//the gap in which the size fits is backwards from here. posFrontofGap returns the position of the end of this contig next to gap.
//does not include N's included in the original contigs, only gaps added by spinner.
{
	ContigInScaffold *nextContig;
	int32my pos=0,nextLen,maxGap,thisGap; //a position holder on scaffold
	//there needs to be a gap in the range, otherwise return
	if (range_min > scaf->len - scaf->backContig->contig->len ||
			range_max <= scaf->frontContig->contig->len ) return(0);
	//range outside the scaffold also breaks the algorithm
	if (range_min < 1) range_min=1 ;
	if (range_max > scaf->len - 1) range_max=scaf->len - 1;
	//advance position pos across contig-gap pairs until in the range 
	nextContig=scaf->frontContig;
	while (pos < range_min) {
		//printf("pos %"PRI32" contig %s gap %"PRI32" ",pos,nextContig->contig->name,nextContig->backwardGap);
		pos += nextContig->contig->len + nextContig->backwardGap; //add on this contig and the following gap.
		nextContig = nextContig->backwards; //move pointer to contig after that gap.
	}
	//we just went past a gap which is now nextContig->forwardGap
	//part of this gap (possibily all of it) is within the range and this part should be counted as a gap in the range.
	thisGap= pos - range_min > nextContig->forwardGap ? nextContig->forwardGap : pos - range_min;
	maxGap=0; //signifies that no gap larger than size has been found.

	//now loop through gap, stopping when the next gap does not start in the range.
	while (pos + nextContig->contig->len < range_max) {
		//printf("pos %"PRI32" contig %s thisgap %"PRI32" ",pos,nextContig->contig->name,thisGap);
		if (thisGap >= size) {
			//printf(" CAND ");
			if (maxGap>0) return(0); //kill off cases where there are two gaps just as good, we don't know which to use.
			maxGap=thisGap;
			*cis=nextContig->forwards; //give the contig from which the relevant gap is backwards
			*posFrontOfGap=pos - nextContig->forwardGap; // and the position at the relevant end of that contig.
		}
		pos += nextContig->contig->len + nextContig->backwardGap; //move over next contig and gap.
		nextContig = nextContig->backwards; //move the contig pointer back one.
		thisGap= nextContig->forwardGap;
	}
	//we just went past a gap which is now nextContig->forwardGap
	//part of this gap (possibly none of it) was still in the range and should be counted.
	thisGap= pos - range_max < 0 ? 0 : thisGap - (pos - range_max);
	//printf("pos %"PRI32" contig %s thisgap %"PRI32" ",pos,nextContig->contig->name,thisGap);
	if (thisGap >= size) {
		if (maxGap>0) return(0); //kill off cases where there are two gaps just as good, we don't know which to use.
		maxGap=thisGap;
		*cis=nextContig->forwards;
		*posFrontOfGap=pos - nextContig->forwardGap;
	}
	return(maxGap >= size);
}

void WriteScaffoldBases(Scaffold *scaf,FILE *fp)
{
	ContigInScaffold *cis;
	char *base,*firstBase;
	int32my i;
	int32my gap;
	for (cis = scaf->frontContig; cis !=NULL; cis=cis->backwards)
	{
		if (cis->dir)
		{
			firstBase=cis->contig->bases;
			for (base=firstBase + cis->contig->len - 1; base >= firstBase; --base)
				fprintf(fp, "%c",ComplementBase(*base));
		}
		else
		{
			fprintf(fp, "%s",cis->contig->bases);
		}
		gap= cis->backwards == NULL ? 0 : ( cis->backwardGap > 0 ? cis->backwardGap : DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP ) ;
		for (i=0;i<gap;++i) fprintf(fp, "N");
	}
	fprintf(fp, "\n");
}

void WriteScaffoldQuals(Scaffold *scaf,FILE *fp)
{
	ContigInScaffold *cis;
	char *qual,*firstQual;
	int32my i;
	int32my gap;
	for (cis = scaf->frontContig; cis !=NULL; cis=cis->backwards)
	{
		if (cis->dir)
		{
			firstQual=cis->contig->quals;
			for (qual=firstQual + cis->contig->len - 1; qual >= firstQual; --qual)
				fprintf(fp, "%c",*qual);
		}
		else
		{
			fprintf(fp, "%s",cis->contig->quals);
		}
	gap= cis->backwards == NULL ? 0 : ( cis->backwardGap > 0 ? cis->backwardGap : DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP ) ;
	for (i=0;i<gap;++i) fprintf(fp, "#");
	}
	fprintf(fp, "\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/* External functions */

void SortScaffoldListByLength(ScaffoldList *scaffoldList)
{
// sorts the linked list of scaffolds so that the largest is the head and so on.
	int32my i;
	Scaffold *scaf, **pscaf;
	Scaffold **scafSortList; 
	
	//get sorted list of remaining scaffolds
	int32my numScafRemaining=0;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next) ++numScafRemaining; //the linked list of remaining scaffolds is counted
	scafSortList = (Scaffold **) MallocWithAssert(sizeof(Scaffold *) * numScafRemaining,"SortScaffoldListByLength, scafSortList");
	for (scaf=scaffoldList->scaffold,i=0; scaf!=NULL; scaf=scaf->next,++i) scafSortList[i]=scaf;
	qsort ( scafSortList, numScafRemaining, sizeof(Scaffold *), CompareScaffoldsByLength ); //sort this list of pointers by the size of the scaffolds
	//set linked list of scaffolds to this order
	scaffoldList->scaffold=*scafSortList; //reset the head of list to the biggest scaffold
	for (pscaf=scafSortList + 1, i=1; i<numScafRemaining -1 ; ++pscaf, ++i) {
		(*pscaf)->next=*(pscaf+1);
		(*pscaf)->prev=*(pscaf-1);
	}
	scafSortList[numScafRemaining - 1]->prev= scafSortList[numScafRemaining - 2];
	scafSortList[numScafRemaining - 1]->next=NULL;
	(*scafSortList)->next=*(scafSortList + 1); //I assume that there are more than 1 scaffold
	(*scafSortList)->prev=NULL;
	free(scafSortList);
}

HalfEdge *GetHalfEdge(Scaffold *scaf1, boolmy_fast dir1, Scaffold *scaf2, boolmy_fast dir2)
{
//  if there is a half edge from scaf1 to scaf2 with given directions, returns a pointer to it, and NULL if not.
//  (as of July 2012 this is only used once, in FindClosest ,and doesn't seem to show much in profiling.)
	HalfEdge *he;
	for (he = dir1 ? scaf1->frontHalfEdge : scaf1->backHalfEdge; he != NULL; he=he->next )
		if (he->oppositeScaffold == scaf2 && he->otherEnd->dir == dir2) break;
	return(he);
}

int32my GetHalfEdgeRange(HalfEdge *he)
{
//  For a given half edge, finds the maximum distance between two good reads supporting the edge
//  (as of July 2012 this plays no role in the scaffolding.  It sometimes get outputted in verbose mode
//  as it might give an indication that some read in the edge is an outlier, or there is something else wrong.
//  at some point it was benig used as a rough estimate of how narrow the distribution of reads for the edge was.)
	Read *read;
	int32my maxpos=0, minpos=2000000000;
	for (read=he->read; read!=NULL; read=read->next)
	{
		if (!read->chimeric && !read->lowQual && !read->crossBiotin) {  //if has not been flagged as suspect
			if (read->scaffoldStart < minpos) minpos=read->scaffoldStart;
			if (read->scaffoldEnd > maxpos) maxpos=read->scaffoldEnd;
		}
	}
	return(maxpos-minpos);
}

void MakeSmallContigsDormant(ScaffoldList *scaffoldList, int smallContigThreshold)
{
// turns on the "dormant" flag for all scaffolds below the threshold size.
// (the intention is that they will not be merged into larger scaffolds at this stage.)
	Scaffold *scaf;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
		if (scaf->len < smallContigThreshold ) scaf->isDormant=1;
}

void ResetDormant(ScaffoldList *scaffoldList)
{
// turns on the "dormant" flag off for all scaffolds.
	Scaffold *scaf;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
		scaf->isDormant=0;
}

float EstimateRelativePosition(Scaffold *scaf, HalfEdge *he)
// does a simple mean estimate of the position of halfEdge->oppositeScaffold relative to the front of scaf
// uses all reads including those marked chimeric and crossbiotin.
// (as of July 2012 used in FindContainedScaffolds in refine.c only).
{
	Read *read1,*read2;
	Scaffold *opScaf;
	float w=0.0, relPos=0.0;
	//float maxWeight=0.0;
	int32my x,endLen1,endLen2, len, opLen;
	//maxWeightStd;
	len=scaf->len;
	opLen=he->oppositeScaffold->len;
	for (read1=he->read;read1!=NULL;read1=read1->next) {
		read2=GetMate(read1);
		//if (read1->crossBiotin) printf("C");
		//printf("R %"PRI32"-%"PRI32" ",read1->scaffoldStart,read2->scaffoldStart); 
		endLen2= read2->dir_s ?  read2->scaffoldStart : (opLen - read2->scaffoldStart) ;
		relPos += (float) read1->lib->weight * (read1->scaffoldStart + (he->dir ?  endLen2 - read1->lib->size : read1->lib->size - endLen2));
		w +=  read1->lib->weight;
		//if (read->lib->weight > maxWeight) {maxWeight=read1->lib->weight;maxWeightStd=read1->lib->std;}
	}
	//printf("\n");
	relPos /= w;
	return (relPos);
}

// for setEdge we need to sort the read pairs supporting the edge by the gap size they imply, in order to do
// things like find quartiles.  That's what these next two things are for.
typedef struct _ReadStat//for sorting read pairs
{
        Read *read1,*read2;
	int32my gap;
} ReadStat;

int32my CompareReadStats(const void * a, const void * b)
{	
	return ( ((ReadStat *)a)->gap - ((ReadStat *)b)->gap );
}

void SetEdge(HalfEdge *halfEdge,Edge *edge,SetEdgeParams *params)
{
//sets the gapSize and all other data in edge, based on read list in halfedge.
	HalfEdge *other;
	Read *read1, *read2;
	int32my i;
	int32my endLen1,endLen2,len1,len2; //to calculate gap sizes
	int32my gap;
	double w,gapSize;
	ReadStat *readArray, *currentRead;
	double q1pos,q1val,q3pos,q3val,iqr,s,s_next,prevGap,nextGap;
	float maxWeight,maxWeightStd;
	PairLib *lib;
	int32my totalScore=0;

	other=halfEdge->otherEnd;
	len1=other->oppositeScaffold->len; //halfedges don't store the scaf they belong to. 
	len2=halfEdge->oppositeScaffold->len; //however they do have the opposite scaffold which we exploit here.
	//allocate an array of the good read pairs to be used for gap estimation.
	readArray  = (ReadStat *) MallocWithAssert(sizeof(ReadStat) * halfEdge->numReads,"SetEdge, readArray");
	//we are going to test the pairs supporting this edge to see which ones are reliable,
	//using three criteria.
	currentRead=readArray;
	w=0.0;
	int numLowQual=0;
	int numChim=0;
	edge->numPair=edge->numGoodPair=edge->numBiotin=0; //counts of various kinds of pairs supporting this edge

	for (read1=halfEdge->read;read1!=NULL;read1=read1->next) //list of reads in matching half edges in same order
	{
		read2=GetMate(read1);

		totalScore += read1->matchScore + read2->matchScore;
		/*(read1->scaffoldStart > len1 *0.333 && read1->scaffoldStart < len1 * 0.667) ? 1 : 0
			+ (read2->scaffoldStart > len1 *0.333 && read2->scaffoldStart < len2 * 0.667) ? 1 : 0;
		*/

		lib=read1->lib;
		//printf("lib size %d lib std %d chim sig %d\n",lib->size,lib->std,params->chimeraSigmas);
		++(edge->numPair);
		endLen1= read1->dir_s ?  read1->scaffoldStart : (len1 - read1->scaffoldStart) ; // gap is pointed to by direction of both reads
		endLen2= read2->dir_s ?  read2->scaffoldStart : (len2 - read2->scaffoldStart) ;
		gap= lib->size - endLen1 - endLen2; //gap implied by this pair assuming the given insert size
		//if it's low quality or inconsistently placed (i.e. the gap size is hugely negative) mark
		//printf("dir %d start %" PRI32 " end %d len %" PRI32 "\n",read1->dir_s,read1->scaffoldStart,len1 - read1->scaffoldStart, len1);
		//printf("dir %d start %" PRI32 " end %d len %" PRI32 "\n",read2->dir_s,read2->scaffoldStart,len2 - read2->scaffoldStart, len2);
		if ( gap < - params->chimeraSigmas * lib->std) //while lowQual and crossBiotin read stusus doesn't change on merges, this does.
			read1->chimeric=read2->chimeric=1;
		else {
			read1->chimeric=read2->chimeric=0;
			if (!read1->lowQual && read1->crossBiotin) ++edge->numBiotin;//num flagged CB but are otherwise okay. Why am I keeping this?
		}
		if (!read1->lowQual && !read1->crossBiotin && !read1->chimeric) {  //if passed checks, update the numGoodPair and add it to our array of read pairs.
			++(edge->numGoodPair);
			currentRead->gap=gap;
			currentRead->read1=read1;
			currentRead->read2=read2;
			++currentRead;
			w += lib->weight; //a number of assumptions here.  May be good to get rid of outliers.
		}
		if (read1->lowQual) ++ numLowQual;
		if (read1->chimeric) ++numChim;
	}
	//edges with enough supporting pairs are labelled using the edgeQual member.
	//everything after that is basically calculating gap size and error.
	//gap size estimation cannot go forward if there are no good pairs,
	// and needn't be done if the edge is invalid because it has not got enough good pairs.
	edge->meanScore= (float) totalScore * 0.5 / (float) halfEdge->numReads;
	edge->edgeQual= edge->numGoodPair >= params->numPairThreshold ? 1 : 0 ;
	//sets a flag that is used in FindClosest, some edges are too weak to join but almost as strong as edges that can join
	//so they should count as conflicting edges. other edges don't need their gap sizes to be estimated.
	if (edge->numGoodPair >= params->numPairThreshold * params->StrengthRatioForWeakEdges && edge->numGoodPair >= 0) edge->conflictCheck=1; else 
	{
		edge->conflictCheck=0;
		edge->w=-1.0;
		edge->gapSize=-1;
		free(readArray);
		return;
	}
	//if there are 4 or less it is not worth taking out outliers using the iqr method.
	if (edge->numGoodPair < 5) { 
		gapSize=0.0;
		w=0.0;
		maxWeight=0.0;
		currentRead=readArray;
		for(i=0; i<edge->numGoodPair; ++i)
		{
			gapSize+= (float) currentRead->gap * currentRead->read1->lib->weight;
			w +=  currentRead->read1->lib->weight;
			if (currentRead->read1->lib->weight > maxWeight) {maxWeight=currentRead->read1->lib->weight;maxWeightStd=currentRead->read1->lib->std;}
			++currentRead;
		}
		gapSize /= w;
		w= w / (maxWeight * maxWeightStd); //rescale w so it can be used as a measure of the uncertainty of the gapSize
		edge->w= w > MIN_INCERTAINTY_REC ? MIN_INCERTAINTY_REC : w;
		edge->gapSize=gapSize;
		//if (gapSize > 150000 || gapSize < -10000 ) printf("(WARNING) SetEdge: w %f gap %f \n",w,gapSize); 
		free(readArray);
		return;
	}
	//otherwise we calculate the quartiles in order to eliminate outliers before taking an average.
	qsort ( readArray, edge->numGoodPair, sizeof(ReadStat), CompareReadStats ); //sort the good reads by the gap size they imply
	//find quartiles (with weights 1/std)
	q1pos= w / 4.0 ;
	q3pos= 3.0 * w / 4.0 ;
	s=0;
	for(i=0; (s_next= s + readArray[i].read1->lib->weight) < q1pos && i<edge->numGoodPair; ++i) //first quartile
		s=s_next;
	if (i==0) q1val=readArray[0].gap; else if (i>=edge->numGoodPair) q1val=readArray[edge->numGoodPair-1].gap; else
		q1val=  readArray[i-1].gap + (float) (readArray[i].gap - readArray[i-1].gap) * (q1pos - s) / readArray[i].read1->lib->weight;

	for( ; (s_next= s + readArray[i].read1->lib->weight) < q3pos && i<edge->numGoodPair; ++i) //carry on to 3rd quartile
		s=s_next;
	if (i==0) q3val=readArray[0].gap; else if (i>=edge->numGoodPair) q3val=readArray[edge->numGoodPair-1].gap; else
		q3val=  readArray[i-1].gap + (float) (readArray[i].gap - readArray[i-1].gap) * (q3pos - s) / readArray[i].read1->lib->weight;
	iqr= q3val - q1val; //interquartile range
	//now go through and estimate gap size eliminating outliers (defined as 1.5 iqr outside of the quartiles).
	gapSize=0.0;
	w=0.0;
	maxWeight=0.0;
	currentRead=readArray;
	//if (iqr<0.1 || iqr > 10000000 ) printf("(WARNING) SetEdge: iqr %f q3val %f q1val %f allocated %d \n",iqr,q3val, q1val,edge->numGoodPair);

	for(i=0; i<edge->numGoodPair; ++i)
	{
		if (currentRead->gap >= q1val - iqr * 1.5) {
			if (currentRead->gap > q3val + iqr * 1.5) break;
			//printf("g %d std %d", currentRead->gap, currentRead->read1->std);
			gapSize+= (float) currentRead->gap * currentRead->read1->lib->weight;
			w +=  currentRead->read1->lib->weight;
			if (currentRead->read1->lib->weight > maxWeight) {maxWeight=currentRead->read1->lib->weight;maxWeightStd=currentRead->read1->lib->std;}
		}
		++currentRead;
	}
	gapSize /= w;
	w= w / (maxWeight * maxWeightStd); //rescale w so it can be used as a measure of the uncertainty of the gapSize
	if (w==0 || ! (gapSize < 10000000) ) printf("(WARNING) SetEdge: w %f gap %f \n",w,gapSize); 
	edge->w=w > MIN_INCERTAINTY_REC ? MIN_INCERTAINTY_REC : w;
	edge->gapSize=gapSize;
	
	//if (gapSize > 150000 || gapSize < -10000 ) printf("(WARNING) SetEdge: w %f gap %f \n",w,gapSize); 
	free(readArray);
}

void QuickEstimateGapNoCutoff(HalfEdge *he)
//finds the gap for the relevant half edgeby weighted average without throwing away "chimeric" edges
//or getting rid of outliers.
{
	Edge *edge;
	Read *read, *readMate;
	float gapSize, w, maxWeight, maxWeightStd;
	int32my gap, len1, len2;

	len1= he->otherEnd->oppositeScaffold->len;
	len2= he->oppositeScaffold->len;
	edge=he->edge;
	gapSize=0.0;
	w=0.0;
	maxWeight=0.0;
	for (read=he->read; read!=NULL; read=read->next) {//list of reads in matching half edges in same order
		if (read->crossBiotin) continue;
		readMate=GetMate(read);
		gap= read->lib->size - (read->dir_s ?  read->scaffoldStart : (len1 - read->scaffoldStart)) -
			(readMate->dir_s ?  readMate->scaffoldStart : (len2 - readMate->scaffoldStart)) ;
		gapSize += (float) gap * read->lib->weight;
		w +=  read->lib->weight;
		if (read->lib->weight > maxWeight) {
			maxWeight= read->lib->weight;
			maxWeightStd=read->lib->std;
		}
	}
	gapSize /= w;
	w= w / (maxWeight * maxWeightStd); //rescale w so it can be used as a measure of the uncertainty of the gapSize
	edge->w= w > MIN_INCERTAINTY_REC ? MIN_INCERTAINTY_REC : w;
	edge->gapSize=gapSize;
}

void RemoveEdge(HalfEdge *he1)
{
//remove the two halfedges (one of which is he1)
// from the lists of edges associated to the two scaffolds
//has to be in properly constructed scaffold list or something odd will happen
	HalfEdge *he2;
	Scaffold *scaf1,*scaf2;
	he2=he1->otherEnd;
	scaf1=he2->oppositeScaffold;
	scaf2=he1->oppositeScaffold;
	//join list over he1
	if (he1->prev == NULL)
		if (he1->dir) scaf1->frontHalfEdge=he1->next; else scaf1->backHalfEdge=he1->next;
	else	he1->prev->next=he1->next;
	if (he1->next != NULL) he1->next->prev=he1->prev; 
	//join list over he2
	if (he2->prev == NULL)
		if (he2->dir) scaf2->frontHalfEdge=he2->next; else scaf2->backHalfEdge=he2->next;
	else	he2->prev->next=he2->next;
	if (he2->next != NULL) he2->next->prev=he2->prev;
	// update number of halfedges in these scaffolds
	if (he1->dir) -- scaf1->frontEdgeNum; else -- scaf1->backEdgeNum;
	if (he2->dir) -- scaf2->frontEdgeNum; else -- scaf2->backEdgeNum;
}

int RemoveScaffold(ScaffoldList *scaffoldList, Scaffold *scaf)
{
//remove the scaffold from the scaffoldlist
//NB: does not remove edges to scaf from active scafs, that should have already happened or be irrelevant for some reason
	if (scaf->prev!=NULL) scaf->prev->next=scaf->next; else scaffoldList->scaffold = scaf->next;
	if (scaf->next!=NULL) scaf->next->prev=scaf->prev;
	return(1);
}

void MakeNameTable(ContigList *contigList)
{
// generate the name table, used to look up the pointers to scaffolds by knowing their name in GetIndexFromName.
// (primarily used to find the correct contig for each read loaded.
//included as part of the loading routines.
	Contig *contig;
	int32my index;

	contigList->nameTable = (NameTable *) MallocWithAssert(sizeof(NameTable)*contigList->contigNum,"MakeNameTable, contigList->nameTable");
	if (contigList->nameTable == NULL) {fprintf(stderr, "LoadFastQ: could not allocate name table\n");exit(1);}
	for (contig=contigList->contig, index=0; index < contigList->contigNum; ++contig, ++index)
	{
		strcpy(contigList->nameTable[index].name,contig->name);
		contigList->nameTable[index].index=index; // good code should be self documenting, like this clear example.
	}
	printf("Sorting for lookup list \n");
	qsort ( contigList->nameTable, contigList->contigNum, sizeof(NameTable), CompareNameTableEntries );
}

int32my GetIndexX(char *name, ContigList *contigList)
{
//returns index to contig of the given name or -1 if not found using nameTable and bsearch.
	NameTable *entry;
	NameTable model;
	strcpy(model.name, name); //this stage won't be the longest anyway, just be lazy
	entry= (NameTable *) bsearch(&model,contigList->nameTable,contigList->contigNum,sizeof(NameTable),CompareNameTableEntries);
	if (entry==NULL) return (-1);
	return(entry->index);
}

int32my GetIndexFromName(char *namex, ContigList *contigList)
{
//returns array index to contig/scaffold of the given name or -1 if not found.
	Contig *contig=NULL;
	if (!strcmp(namex,"*")) return(-1); //non-aligned reads
	if (contigList->numInName > -1) return (contigList->numToIndex[atoi(namex+contigList->numInName)] );
	else return(GetIndexX(namex, contigList));
}

Contig *GetContigFromName(char *namex, ContigList *contigList)
{
//returns pointer to contig of the given name.
	int32my index;
	index=GetIndexFromName(namex, contigList);
	if (index == -1)  return (NULL);
	return(contigList->contig + index);
}

void SetScaffoldIndex(Read *read, char *contigName, ContigList *contigList)
{
// sets read->scaffoldIndex to the index of the scaffold on which read 1 is matched and sets the variable in the read stucture to the same value.
//returns -1 if the name is not in the list (this happens sometimes as some contigs are not good and are not included in fastq file).
//  No need for names to be in any order but this weorks much faster if numInName is supplied.
//  this only needs to be called before scaffolding commences.
	read->scaffoldIndex=GetIndexFromName(contigName, contigList);
}

void CrossBiotinFlag(HalfEdge *he, int16f SWdeficit)
{
//checks if reads supporting the half edge are crossBiotin and flags them accordingly
	Read *read1,*read2;
	//go through reads and the paired read on the matching halfedge.
	for (read1=he->read;read1!=NULL;read1=read1->next) {
		read2=GetMate(read1);
		read1->crossBiotin=read2->crossBiotin=IsCrossBiotin(read1,read2,SWdeficit);
	}
}

void LowMatchScoreFlag(HalfEdge *he, char threshold)
{
//checks if reads supporting the half edge are low match score and flags them accordingly
	Read *read1,*read2;
	//go through reads and the paired read on the matching halfedge.
	for (read1=he->read;read1!=NULL;read1=read1->next) {
		read2=GetMate(read1);
		read1->lowQual = read2->lowQual = (read1->matchScore <= threshold || read2->matchScore <= threshold);
	}
}

boolmy_fast LoadFastq(char *fname, ContigList *contigList)
{
// fills in config info from fastQ file: names, lengths and base and qual info
// allocates memory block pointed to by contigList->fastQData
// NB: having two different routines for fasta and fastq here is obviously a design mistake.
// But hopefully these won't need much alteration.
// here also used to be a version that loaded everything but the base and qual strings,
// which was useful for estimating inserts.  It is gone now.

	int64my length, i, nBlock;
	int32my blockSize, j, remain;
	char *FQdata;
	char *ptr, *ptr2;
	FILE *fp;
	int32my numContig, fnLen;
	char *contigName;
	int32my contigLen, contigLen2;
	Contig *contig, *hash=NULL;
	int32my k;
	int32my mult, num;
	int32my *iptr, numToIndexSize;
	
	fp = fopen(fname,"r");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}
	fseek(fp, 0, SEEK_END);
	length = ftell(fp);
	rewind(fp);

	blockSize = 1024;
	FQdata = (char *) MallocWithAssert(sizeof(char)*length,"LoadFastq, FQdata");
	nBlock = length/blockSize;
	remain = length - nBlock*blockSize;
	if((j=fread(FQdata, blockSize, nBlock, fp)) != nBlock){
		fprintf(stderr,"file reading %s error: block reading\n",fname);
		exit(1);
	}
	if(remain != 0){
		if((j=fread(FQdata+blockSize*nBlock, sizeof(char), remain, fp)) != remain){
			fprintf(stderr,"file reading %s error: remain reading\n",fname);
			exit(1);
		}
	}

        numContig = 0;
        ptr = FQdata;
        for(i=0;i<length;++i) {
                if(*ptr == '@' && (i==0 || *(ptr-1) == '\n') ) ++numContig; //counts all subsequent contigs
		++ptr;
	}
	contigList->contigNum = numContig;
	printf("----- LoadFastq: Contains data for %d contigs \n",contigList->contigNum);

        //---allocate mem for contigList------------------------
        if (contigList->contig!=NULL) {fprintf(stderr,"LoadContigListFromFastq: contiglist->contig not NULL \n"); exit(1);}
        contigList->contig  = (Contig *) MallocWithAssert(sizeof(Contig) * numContig,"LoadContigListFromFastq, contigList->contig");
	if (contigList->numInName > -1) {
		numToIndexSize=numContig * 2;
		contigList->numToIndex = (int32my *)  MallocWithAssert(sizeof(int32my) * numToIndexSize,"LoadContigListFromFastq, contigList->NumToIndex");
		for (i=0, iptr=contigList->numToIndex; i < numToIndexSize; ++i, ++iptr) *iptr=-1;
	}
        //printf("done allocation \n");
        contig=contigList->contig;
	ptr = FQdata;
	k = 0;
	while(ptr-FQdata < length){
		if(*ptr++ != '@'){
			fprintf(stderr,"LoadFastq: No \'@\' when expected, number of contigs read: %d \n",k);
			exit(1);
		}
		contigName = ptr;
		fnLen = 0;
		while(*ptr++ != '\n'){
			++fnLen;
		}
		ptr[-1] = '\0';
		if (fnLen > ReadNameLen) {fprintf(stderr,"LoadContigListFromFastq: contig name too long \n"); exit(1);}
		strcpy(contig->name, contigName);
		if (contigList->numInName > -1) {
			num=atoi(contigName + contigList->numInName );
			if (num >= numToIndexSize) {
				if (num > 10000000) {fprintf(stderr,"LoadContigListFromFastq: contig number from name is too big:  %d \n",num); exit(1);}
				mult=0;
				while((numToIndexSize * mult) <= num) ++mult;
				contigList->numToIndex=(int32my *) ReallocWithAssert(contigList->numToIndex, sizeof(int32my) * numToIndexSize*mult,"LoadContigListFromFastq, contigList->contig");
				for (i=numToIndexSize, iptr=contigList->numToIndex + numToIndexSize ; i < numToIndexSize*mult; ++i, ++iptr) *iptr=-1;
				numToIndexSize=numToIndexSize*mult;
			}
			contigList->numToIndex[num]=k;
		}
		contig->bases = (char *)ptr;
		ptr2 = ptr;
		contigLen=0;
		while(*ptr != '+'){
			if(*ptr == '\n') {
				ptr++;
				continue;
			}
			*ptr2++ = *ptr++;
			contigLen++;
		}
		*ptr2 = '\0';
		contig->len = contigLen;

		while(*ptr != '\n') ++ptr; //sometimes the name is above the quals, sometimes not, best just to advance to the newline
		if(*ptr != '\n'){
			fprintf(stderr,"FQread error-4 \n");
			exit(1);
		}
		*ptr++ = '\0';
		contig->quals = (char *)ptr;
		contigLen2 = 0;
		ptr2 = ptr;
                while( ptr != FQdata+length && (contigLen2 < contigLen || *ptr != '@') ){
			//
                        if(*ptr == '\n') {
                                ptr++;
                                continue;
                        }
                        *ptr2++ = *ptr++; //this kind of stuff was there in RPono.
			contigLen2++;
                }
		*ptr2 = '\0';
		if(contigLen2 != contigLen){
	                fprintf(stderr,"Warning: different length for base and qual: %ld %ld\n", contigLen,  contigLen2);
		}
		++k;
		++contig;
	}
	printf("----- LoadFastq: Read fastq file.\n",k);
	contigList->FastQData=FQdata;
	contigList->FastQDataRetained=1;
	//for (i=0;i<contigList->contigNum;++i)
	//	fprintf(stderr,"contigList->numToIndex[%"PRI64"]=%"PRI32"\n",i,contigList->numToIndex[i]);

	if (contigList->numInName < 0) MakeNameTable(contigList);
	return (1);
}

boolmy_fast LoadFasta(char *fname, ContigList *contigList)
{
// fills in config info from fasta file: names, lengths and base info (no quals in this case)
// allocates memory block pointed to by contigList->fastQData
// note:

	int64my length, i, nBlock;
	int32my blockSize, j, remain;
	char *FQdata;
	char *ptr, *ptr2;
	FILE *fp;
	int32my numContig, fnLen;
	char *fname1, *fname2;
	int32my contigLen, contigLen2;
	Contig *contig, *hash;
	int32my k;
	int32my mult, num;
	int32my *iptr, numToIndexSize;
	
	fp = fopen(fname,"r");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}
	fseek(fp, 0, SEEK_END);
	length = ftell(fp); //+1 in case the file ends without a \n, I need to add a \0 in the final slot.
	rewind(fp);

	blockSize = 1024;

	FQdata = (char *) MallocWithAssert(sizeof(char)*length + 1, "LoadFasta, FQdata");
	nBlock = length/blockSize;
	remain = length - nBlock*blockSize;

	if((j=fread(FQdata, blockSize, nBlock, fp)) != nBlock){
		fprintf(stderr,"file reading %s error: block reading\n",fname);
		exit(1);
	}

	if(remain != 0){
		if((j=fread(FQdata+blockSize*nBlock, sizeof(char), remain, fp)) != remain ){
			fprintf(stderr,"file reading %s error: remain reading\n",fname);
			exit(1);
		}
	}

	numContig = 0;
        ptr = FQdata;
        for(i=0;i<length;++i) {
                if(*ptr == '>' && (*(ptr-1) == '\n' || i==0) ) ++numContig; //counts all contigs
		++ptr;
	}
	contigList->contigNum = numContig; 
	printf("Contains data for %d contigs \n",contigList->contigNum);

        //---allocate mem for contigList------------------------
        if (contigList->contig!=NULL) {fprintf(stderr,"LoadFasta: contiglist->contig not NULL \n"); exit(1);}
        contigList->contig  = (Contig *) MallocWithAssert(sizeof(Contig) * numContig,"LoadFasta, contigList->contig");
	if (contigList->numInName > -1) {
		contigList->numToIndex = (int *)  MallocWithAssert(sizeof(int) * numContig * 2,"LoadContigListFromFastq, contigList->NumToIndex");
		for (i=0, iptr=contigList->numToIndex; i<numContig * 2; ++i, ++iptr) *iptr=-1;
		numToIndexSize=numContig * 2;
	}
        //printf("done allocation \n");
        contig=contigList->contig;
	ptr = FQdata;
	k = 0;
	if(*ptr != '>'){
			fprintf(stderr,"LoadFasta: - no > at beginning of file \n");
			exit(1);
	}
	while(ptr-FQdata < length){
		++ptr; //gets over the '>'
		fname1 = ptr;
		fnLen = 0;
		while(*ptr++ != '\n'){
			++fnLen;
		}
		ptr[-1] = '\0';
		if (fnLen > ReadNameLen) {fprintf(stderr,"LoadFasta: contig name too long \n"); exit(1);}
		strcpy(contig->name, fname1);
		if (contigList->numInName > -1) {
			num=atoi(fname1 + contigList->numInName );
			if (num >= numToIndexSize) {
				if (num > 100000000) {fprintf(stderr,"LoadContigListFromFastq: contig number from name is too big:  %d \n",num); exit(1);}
				mult=0;
				while((numToIndexSize * mult) < num) ++mult;
				contigList->numToIndex=ReallocWithAssert(contigList->numToIndex, sizeof(int) * numToIndexSize * mult,"LoadContigListFromFastq, contigList->contig");
				for (i=numToIndexSize, iptr=contigList->numToIndex + numToIndexSize; i < numToIndexSize * mult; ++i, ++iptr) *iptr=-1;
				numToIndexSize=numToIndexSize*mult;
			}
			contigList->numToIndex[num]=k;
			//printf("%s %"PRI32" %"PRI32" ",fname1,num,k);
		}
		//printf("name of contig %d is ~%s~ \n",k,contig[k].name);
		contig->bases = ptr;
		ptr2 = ptr;
		contigLen=0;
		while(ptr != FQdata+length && *ptr != '>') //I suppose this strips line breaks that might sometimes occur before the next entry
			if(*ptr == '\n')
				ptr++;
			else {
				*ptr2++ = *ptr++;
				contigLen++;
			}
		*ptr2 = '\0';
		contig->len = contigLen;
		contig->numReads=0;
		++k;
		++contig;
	}
	printf("Read data for %d contigs \n",k);
	contigList->FastQData=FQdata;
	contigList->FastQDataRetained=1;
	if (contigList->numInName < 0) MakeNameTable(contigList);
	return (1);
}

int SaveFastq(ScaffoldList *scaffoldList, char *fname)
{
	FILE *fp;
	Scaffold *scaf;
	int32my i;

	printf("Saving in fastq format\n");
	fp = fopen(fname,"w");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}

	for (scaf=scaffoldList->scaffold,i=0;scaf!=NULL;scaf=scaf->next,++i)
	{
		fprintf(fp, "@%s\n", scaf->frontContig->contig->name);  //original name of the front contig in scaf will be name of scaffold
		WriteScaffoldBases(scaf,fp);
		fprintf(fp, "+%s\n", scaf->frontContig->contig->name);
		WriteScaffoldQuals(scaf,fp);
	}
	printf("Saved %d contigs\n",i);
	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
}

int SaveFasta(ScaffoldList *scaffoldList, char *fname)
{
	FILE *fp;
	Scaffold *scaf;
	int32my i;

	printf("Saving in fasta format\n");
	fp = fopen(fname,"w");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}
	for (scaf=scaffoldList->scaffold,i=0;scaf!=NULL;scaf=scaf->next,++i)  {
		fprintf(fp, ">%s\n", scaf->frontContig->contig->name);  //original name of the front contig in scaf will be name of scaffold
		WriteScaffoldBases(scaf,fp);
	}
	printf("Saved %d contigs\n",i);
	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
}

int SaveContigsInScaffoldsProperFormat(ScaffoldList *scaffoldList, char *fname, int dir)
{
// saves a file with information about what contigs are in what scaffolds, etc.
	FILE *fp;
	Scaffold *scaf;
	ContigInScaffold *cis;
	int64my lenSum, gapSum; //the cumulative length up to the start of that gap (useful when finding bad joins)
	int32my gap;
	char mark;

	printf("Saving contigs in file %s\n",fname);
	fp = fopen(fname,"w");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
	{
		cis= dir ? scaf->frontContig : scaf->backContig;
		fprintf(fp, "supercontig %s len %" PRI32 "junction  %d<-->%d\n",
				 cis->contig->name, scaf->len, scaf->isFrontJunction, scaf->isBackJunction );
		fprintf(fp, "contig %s len %" PRI32 " reads/bp %f\n", 
				cis->contig->name, cis->contig->len, 
				cis->contig->numReads / ((float) (cis->contig->len - 100 > 0 ? cis->contig->len - 100 : 1)));
		lenSum=cis->contig->len;
		gapSum=0;
		for (cis=dir ? cis->backwards: cis->forwards ; cis!=NULL; cis=dir ? cis->backwards: cis->forwards) {
			gap=dir ? cis->forwardGap : cis->backwardGap;
			mark=dir ? cis->forwardMark : cis->backwardMark;
			fprintf(fp, "gap %" PRI32 " %d %" PRI64 " %d\n",
				gap,
				(int) dir ? cis->forwardHe->edge->numGoodPair : cis->backwardHe->edge->numGoodPair,
				lenSum,
				(int) mark );
			lenSum += cis->contig->len + (gap > 0 ? gap : DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP);
			gapSum += gap > 0 ? gap : DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP;
			fprintf(fp, "contig %s len %" PRI32 " reads/bp %f\n", cis->contig->name, cis->contig->len, (float) cis->contig->numReads / ((float) (cis->contig->len - 100 > 0 ? cis->contig->len - 100 : 1) ));
		}
		fprintf(fp, "\n");
	}
	printf("(GAP) total length of added gaps: %ld \n", gapSum);
	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
}

void ClearMarks(ScaffoldList *scaffoldList)
//Mark the first n connected components containing the scafs in order
{
	Scaffold *scaf;
	for (scaf=scaffoldList->scaffold; scaf!=NULL; scaf=scaf->next)
		scaf->mark=0;
}

void MarkConnectedComponent(Scaffold *scaf, int mark, int depth)
{
//part of SaveGDF, a way of reducing number of nodes saved.
// this doesn't really give all nodes of the given depth away from the start point as the search is depth first.
	HalfEdge *he;
	if (scaf->mark || ! depth ) return;
	scaf->mark=mark;
	for (he= scaf->frontHalfEdge; he !=NULL ; he = he->next)
		if (he->edge->numGoodPair >= 5) MarkConnectedComponent(he->oppositeScaffold, mark, depth -1);
	for (he= scaf->backHalfEdge; he !=NULL ; he = he->next)
		if (he->edge->numGoodPair >= 5) MarkConnectedComponent(he->oppositeScaffold, mark, depth -1);
}

void MarkConnectedComponents(ScaffoldList *scaffoldList, int num, int depth)
//Mark the first n connected components containing the scafs in order
{
	Scaffold *scaf;
	int32my mark=0;

	for (scaf=scaffoldList->scaffold; scaf!=NULL && mark<num ; scaf=scaf->next) {
		if (scaf->mark==0) {
			++mark;
			MarkConnectedComponent(scaf, mark, depth);
		}
	}
}

int SaveGDF(ScaffoldList *scaffoldList, char *fname, int min_size, int num, int depth)
{
// saves a graph file showing the scaffolding graph
//if num==0 saves all, if > 0 saves num connected components only
	FILE *fp;
	Scaffold *scaf;
	HalfEdge *he;
	int64my pos;
	int64my totalGap;
	char *color, red[]="\"255,0,0\"", green[]="\"0,255,0\"",black[]="\"0,0,0\"";
	Read *read,*read2;

	fp = fopen(fname,"w");
	if (fp == NULL) {fprintf(stderr, "cannot open file %s\n", fname);exit(1);}
	//if requested find connected components
	if (num > 0) {
		ClearMarks(scaffoldList);
		MarkConnectedComponents(scaffoldList, num, depth);
	}
	//write node info
	fprintf(fp, "nodedef>name VARCHAR,label VARCHAR, len INT\n");
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next) {
		if (scaf->len < min_size) continue;
		if (num > 0 && scaf->mark==0) continue;
		fprintf(fp, "%s,%s,%d\n", scaf->frontContig->contig->name, scaf->frontContig->contig->name, scaf->len);
	}
	//write edge info
	fprintf(fp, "edgedef>node1 VARCHAR,node2 VARCHAR, color VARCHAR, gap DOUBLE, goodpairs INT\n");
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next) {
		if (scaf->len < min_size) continue;
		if (num > 0 && scaf->mark==0) continue;
		for (he=scaf->backHalfEdge; he!=NULL ; he=he->next) {
			//if (! he->edge->edgeQual) continue;
			if (he->oppositeScaffold->len < min_size) continue;
			if ( he->edge->numGoodPair < 5 || (num > 0 && ! he->oppositeScaffold->mark) ) continue;
			color = he->edge->edgeQual ? green : black ; //(he == scaf->frontClosest || he == scaf->backClosest)
			fprintf(fp, "%s,%s,%s,%f,%d\n", scaf->frontContig->contig->name,
				he->oppositeScaffold->frontContig->contig->name,color, he->edge->gapSize, he->edge->numGoodPair);
		}
		for (he=scaf->frontHalfEdge; he!=NULL ; he=he->next) {
			//if (! he->edge->edgeQual) continue;
			if (he->oppositeScaffold->len < min_size) continue;
			if ( he->edge->numGoodPair < 5 || (num > 0 && ! he->oppositeScaffold->mark) ) continue;
			color = he->edge->edgeQual ? green :  red ; //(he == scaf->frontClosest || he == scaf->backClosest)
			fprintf(fp, "%s,%s,%s,%f,%d\n", scaf->frontContig->contig->name,
				he->oppositeScaffold->frontContig->contig->name,color, he->edge->gapSize, he->edge->numGoodPair);
		}
	}
	if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);};
}

void InitialiseContigList(ContigList *contigList)
{
// sets num of entries in list to zero and makes list empty
	contigList->contigNum=0;
	contigList->contig=NULL;
	contigList->FastQDataRetained=0;
}

ScaffoldList *SetScaffoldListFromContigList(ContigList *contigList,ScaffoldList *scaffoldList)
//intialises scaffoldList, creating one scaffold for each contig containing only that contig
// there seems to be a lot of redeundant intitalisation here, I don't know why I did that.
{
	ContigInScaffold *CSArray; //the initial locations of the members of the lists of contigs making up scaffolds.  
	ContigInScaffold *contigInScaffold;
	Contig *contig;
	Scaffold *scaffold;
	int32my num; // number of contigs (= number of initial scaffolds)
	int32my i;

	scaffoldList->firstHalfEdge=NULL;
	scaffoldList->firstEdge=NULL; 

	num=contigList->contigNum; 
	scaffoldList->scaffoldNum=num;// initially there are as many scaffolds as contigs, and there will only ever be less (or so I plan)

	scaffoldList->scaffold  = scaffoldList->firstScaffold = //so allocate enough memory for scaffolds
		(Scaffold *) MallocWithAssert(sizeof(Scaffold) * num, "SetScaffoldListFromContigList, scaffoldList->scaffold");  
	//at this stage I will allow more contigs in scaffolds than there are contigs in case I later want to account for repeats etc.
	CSArray=(ContigInScaffold *) MallocWithAssert(sizeof(ContigInScaffold) * num * 2, "SetScaffoldListFromContigList, CSArray");
	scaffoldList->firstContigInScaffold=CSArray; //helps with deallocation later on.
	scaffoldList->newContigInScaffold=CSArray+num;  //where to add new ones later on.
	// initially the scaffolds will be stored in the same order as contigs, although, unlike the contig array, scaffoldlist is later to be used as a linked list.
	
	//fill in the contig in scaffold and contig info
	for (i=0;i<num;++i) {
		contigInScaffold=CSArray+i;	//get the structures for this iteration of loop
		contig=contigList->contig+i;			
		scaffold=scaffoldList->scaffold+i;

		//first fill in contig in scaffold
		//not sure how things are intialised by compiler but I assume they are full of junk.
		contigInScaffold->contig=contig;
		contigInScaffold->dir=0;
        	contigInScaffold->forwards=NULL;    //this single contig is front and back of list intially
        	contigInScaffold->backwards=NULL;
		contigInScaffold->forwardGap=0;
		contigInScaffold->backwardGap=0;
		//now do the scaffold	
		scaffold->mark=0;
		scaffold->frontEdgeNum=0;
		scaffold->backEdgeNum=0;
		scaffold->frontHalfEdge=NULL;
		scaffold->backHalfEdge=NULL; //no graph info yet (not sure how these things get intitialised unless I do it explicitly)		
		scaffold->len=contig->len;                          //length of scaffold = length of contig initially
		scaffold->numReads=contig->numReads;
		scaffold->contigNum=1;                    //number of contigs in this scaffold
		scaffold->frontContig=CSArray+i;   //there's one contiginscaffold in the list initially and they are simply ordered at this point->
		scaffold->backContig=CSArray+i;
		scaffold->isDormant=0;
		scaffold->next=scaffold+1;           //points to next scaffold in list of scaffolds(order not important)
		scaffold->prev=scaffold-1;
	}
	scaffold->next=NULL; //terminate list
	scaffoldList->scaffold->prev=NULL;
	return(scaffoldList);
}

////////// the following are for MakeEdgesFromReads, the routine that take the reads and works out the edges and their properties from them.

typedef struct _ReadComp//for sorting read pairs
{
	Read *read;
	int32my opScafIndex;
	boolmy_fast dir,opDir;
} ReadComp;

int CompareReadComp(const void *a, const void *b)
//half edges are INITIALLY ordered first by the opposite scaffold and then by the dir of the opposite half edge
//if the half edge is larger than the supplied data then +1, smaller -1, equal 0. 
{
	if (((ReadComp *)a)->opScafIndex < ((ReadComp *)b)->opScafIndex) return(-1);
	if (((ReadComp *)a)->opScafIndex > ((ReadComp *)b)->opScafIndex) return(1);
	if (((ReadComp *)a)->dir < ((ReadComp *)b)->dir) return(-1);
	if (((ReadComp *)a)->dir > ((ReadComp *)b)->dir) return(1);
	if (((ReadComp *)a)->opDir < ((ReadComp *)b)->opDir) return(-1);
	if (((ReadComp *)a)->opDir > ((ReadComp *)b)->opDir) return(1);
	//if (((ReadComp *)a)->read < ((ReadComp *)b)->read) return(-1);
	//if (((ReadComp *)a)->read > ((ReadComp *)b)->read) return(1);
	return(0);

}

int CompareReadCompWithoutReadCheck(const void *a, const void *b)
//half edges are INITIALLY ordered first by the opposite scaffold and then by the dir of the opposite half edge
//if the half edge is larger than the supplied data then +1, smaller -1, equal 0. 
{
	if (((ReadComp *)a)->opScafIndex < ((ReadComp *)b)->opScafIndex) return(-1);
	if (((ReadComp *)a)->opScafIndex > ((ReadComp *)b)->opScafIndex) return(1);
	if (((ReadComp *)a)->dir < ((ReadComp *)b)->dir) return(-1);
	if (((ReadComp *)a)->dir > ((ReadComp *)b)->dir) return(1);
	if (((ReadComp *)a)->opDir < ((ReadComp *)b)->opDir) return(-1);
	if (((ReadComp *)a)->opDir > ((ReadComp *)b)->opDir) return(1);
	return(0);

}

int32my SortReads(ReadList *readList, Read **firstRead, int32my numReads, ReadComp *comp)
{
//sorts the reads and also returns number of halfedges needed for this scaffold.
	Read *read,*opRead;
	ReadComp *cp;
	int32my i,numHe;
	if (*firstRead==NULL) return(0);
	for (read=*firstRead, cp=comp, i=0; i < numReads; read=read->next, ++cp, ++i) {
		cp->read=read;
		cp->opScafIndex=(opRead = GetMate(read))->scaffoldIndex;
		cp->dir=read->dir_s;
		cp->opDir=opRead->dir_s;
	}
	//if (read!=NULL) {fprintf(stderr, "pob \n"); exit(1);}
	qsort ( comp, numReads, sizeof(ReadComp), CompareReadComp );
	numHe=1; //we know there is at least one read in the list so there must be one.
	for (i=1,cp=comp+1;i < numReads;++cp, ++i)
		if(CompareReadCompWithoutReadCheck(cp,cp-1)) ++numHe; //if this read will belong to a different halfedge from the previous one.
	*firstRead=comp->read; //reset read list head
	for (cp=comp, i=0; i < numReads - 1; ++cp, ++i)
		cp->read->next=(cp + 1)->read;
	cp->read->next=NULL;
	return(numHe);
}

HalfEdge *AddHalfEdgeToList(HalfEdge *he, Scaffold *scaf, boolmy_fast dir)
{
// adds a halfedge (pre-allocated) to the doubly linked list of he's (for the given scaffold)

	if (dir) {
		if (scaf->frontHalfEdge == NULL) { //if the list of half edges is empty
			scaf->frontHalfEdge=he;
			he->prev=he->next=NULL;
		} else {
			scaf->frontHalfEdge->prev=he;
			he->next= scaf->frontHalfEdge;
			scaf->frontHalfEdge=he;
			he->prev=NULL; //I suppose this is already true but anyway.
		}
		++ scaf->frontEdgeNum;
	} else { //same in the back list
		if (scaf->backHalfEdge == NULL) {
			scaf->backHalfEdge=he;
			he->prev=he->next=NULL;
		} else {
			scaf->backHalfEdge->prev=he;
			he->next= scaf->backHalfEdge;
			scaf->backHalfEdge=he;
			he->prev=NULL;
		}
		++ scaf->backEdgeNum;
	}
}

void SetHalfEdgeData(HalfEdge *he, Scaffold *opScaf, HalfEdge *otherEnd, boolmy_fast dir, Read *read, int32my numReads)
{
//sets the properties of the half edge.
	he->oppositeScaffold=opScaf;
	he->otherEnd=otherEnd;
	he->dir=dir;
	he->read=read;
	he->numReads=numReads;
}

ScaffoldList *MakeEdgesFromReads(ReadList *readList, ContigList *contigList,ScaffoldList *scaffoldList)
{
// uses the readlist to complete the intial data structures:
// the "edges" (structures containing info about the read pairs that connect two scaffolds, when there are any)
// are allocated, their properties set (including a list of supporting reads), and they are listed
// in the scaffolds to which they are associated.
// The ScaffoldList should already have been intialised (with each scaffold containing one contig and no edges).
// a bidirectional edge is represented by 3 structures: 2 half edges and an edge.
// (note: this sorting-based implementation replaces a shorter piece of code that 
// searched all halfedges for each read.  This version runs much more quickly
// when there is a large set of reads).
	Read *read, *nextRead, *firstRead, *lastRead, *opRead, *x;
	HalfEdge *he;
	Edge *edge;
	Scaffold *scaf, *opScaf;
	int32my opScafIndex;
	boolmy_fast dir,opDir;	
	int64my numHe, i;
	int32my numReads, scafIndex, maxNumReads;
	int32my *readNum, *p;
	ReadComp *comp;
	Read **scafRead, **listHead;
	Read oldHead;
	int64my numHeFound,numHeCreated, numBadName;

	readNum= (int32my *) MallocWithAssert(sizeof(int32my) * scaffoldList->scaffoldNum, "MakeEdgesFromReads, readNum" ); 
	scafRead= (Read **) MallocWithAssert(sizeof(Read *) * scaffoldList->scaffoldNum, "MakeEdgesFromReads, nscafRead" );
	for (i=0, p=readNum, listHead=scafRead;i<scaffoldList->scaffoldNum;++i,++p,++listHead) {*p=0; *listHead=NULL;}
	//build list of all reads associated to each scaf (NOT the use of the read->next pointer outside this function)
	for(i=0,read=readList->read; i<readList->readNum; ++i,++read) {
		if ( (scafIndex=read->scaffoldIndex)!=-1 && GetMate(read)->scaffoldIndex!=-1 ) {
			if ( scafIndex >= scaffoldList->scaffoldNum ) {fprintf(stderr, "prob\n"); exit(1);}
			read->next=*(listHead=scafRead + scafIndex);
			*listHead=read;
			++readNum[scafIndex];
		}
	}
	maxNumReads=0;
	for (i=0, p=readNum;i<scaffoldList->scaffoldNum;++i,++p)
		if (*p > maxNumReads) maxNumReads=*p;
	comp = (ReadComp *) MallocWithAssert(sizeof(ReadComp)*maxNumReads, "MakeEdgesFromReads, newHalfEdge" ); 
	numHe=0;
	for (i=0, p=readNum, listHead=scafRead;i<scaffoldList->scaffoldNum;++i,++p,++listHead)
		numHe+=SortReads(readList, listHead, *p, comp); //sorts according to scaffold, opposite scaffold, direction and opposite direction.
	free(comp);
	printf("# Half edges: %d \n", numHe);
	he = (HalfEdge *) MallocWithAssert(sizeof(HalfEdge)*numHe, "MakeEdgesFromReads, he" ); 
	scaffoldList->firstHalfEdge=he; //for deallocation later
	edge = (Edge *) MallocWithAssert(sizeof(Edge)*(numHe/2), "MakeEdgesFromReads, edge" );
	scaffoldList->firstEdge=edge;
	scaffoldList->edgeNum=numHe/2;

	numBadName=numHeFound=numHeCreated=0;
	for (i=0, listHead=scafRead;i<scaffoldList->scaffoldNum;++i,++listHead) {
		read=*listHead;
		while(read!=NULL) {
			firstRead=read;
			opScafIndex=(opRead=GetMate(read))->scaffoldIndex;
			dir=read->dir_s;
			opDir=opRead->dir_s;
			numReads=0;
			while ( read!=NULL && opScafIndex==(x=GetMate(read))->scaffoldIndex &&  dir==read->dir_s && opDir==x->dir_s) {
				++numReads;
				lastRead=read;
				read=read->next;
			}
			lastRead->next=NULL; //break the scaffold's read list into lists for each half edge.
			++numHeFound;
			if (opScafIndex==-1 || opScafIndex >= scaffoldList->scaffoldNum ) ++numBadName;
			if (opScafIndex > i) {//create half edges in this case, already created (or it's -1) otherwise (no == case should be found)
				scaf=scaffoldList->scaffold + i; //i is the index of this scaffold and scaffoldList->scaffold is the array
				opScaf=scaffoldList->scaffold + opScafIndex;
				AddHalfEdgeToList(he,scaf,dir);
				SetHalfEdgeData(he, opScaf, he+1, dir, firstRead, numReads);  //half edge stores these. he+1 is the opposite half edge.
				++he;
				AddHalfEdgeToList(he,opScaf,opDir);
				SetHalfEdgeData(he, scaf, he-1, opDir, opRead, numReads);
				++he;
				numHeCreated+=2;
			}
		}
	}
	printf("allocated for %" PRI64 " half edges, found %" PRI64 " and created %" PRI64 ", with %" PRI64 " bad name \n",
		numHe, numHeFound,numHeCreated,numBadName);

	printf("allocated %" PRI64 " bytes for %" PRI64 " edges (would have been %" PRI64 " without checking number) \n",
		(int64my) sizeof(Edge)*numHe/2,numHe/2, (int64my) sizeof(Edge)*readList->readNum);
	printf("allocated %" PRI64 " bytes for %" PRI64 " half edges (would have been %" PRI64 " without checking number) \n",
		(int64my) sizeof(HalfEdge)*numHe, numHe, (int64my) sizeof(HalfEdge)*readList->readNum);

	//copies of the same read pairs are sometimes present
	//these would confuse the counting of read pairs confirming an edge, and gap estimates, unless removed
	printf("Removing duplicate pairs \n"),
	RemoveDuplicatePairs(scaffoldList);

	//now that we have the edges in memory it's time to fill in their data.
	//this data depends on the reads again.  Now the relevant reads (- duplicates)
	//are stored in lists in the half edges and so are easy to find.
	for (i=0;i<numHe; i+=2)  //halfedges are next to their other ends in the list at this stage.  Only need one of pair.
		AddEdge(scaffoldList->firstHalfEdge+i,&edge); // makes edge, adds pointers from both half edges.

	printf("Number of edges: %d \n",scaffoldList->edgeNum);
	//now the graph in scaffoldList,including all half edges and edges, is ready to go.
	free(readNum);
	free(scafRead);
	return(scaffoldList);
}



void DisplayClosestScaffolds(ScaffoldList *scaffoldList)
{
//prints out a list of scaffolds together with their front and back closest and directional properties
	Scaffold *scaf;
	char *name1,*name2;
	boolmy_fast bdir1,bdir2,fdir1,fdir2;
	char junction[]="junction------junction\0";
	char no_link[]="no_link--------no_link\0";
	

	printf("(CLOSEST) Displaying estimated closest scaffolds \n");;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
	{
		name1=scaf->frontClosest != NULL ? scaf->frontClosest->oppositeScaffold->frontContig->contig->name : 
			(scaf->isFrontJunction ? junction : no_link ) ;
		name2=scaf->backClosest != NULL ? scaf->backClosest->oppositeScaffold->frontContig->contig->name : 
			(scaf->isBackJunction ? junction : no_link ) ;
		fdir1=scaf->frontClosest != NULL ? scaf->frontClosest->dir : 9;
		bdir1=scaf->backClosest != NULL ? scaf->backClosest->dir : 9;
		fdir2=scaf->frontClosest != NULL ? scaf->frontClosest->otherEnd->dir : 9;
		bdir2=scaf->backClosest != NULL ? scaf->backClosest->otherEnd->dir : 9;
		printf(" %s <%d--%d- fC %s bC %s -%d--%d> %s len: %d dormant? %d\n", 
			name1,fdir2,fdir1,scaf->frontContig->contig->name,scaf->backContig->contig->name,bdir1,bdir2,name2,scaf->len,scaf->isDormant);
		//printf(" %d <- %d -> %d \n", len1,scaf->len,len2);
	}
}

int DeleteScaffoldList(ScaffoldList *scaffoldList)
{
//free memory in scaffoldList
	free(scaffoldList->firstScaffold);
	free(scaffoldList->firstHalfEdge);
	free(scaffoldList->firstEdge);
	free(scaffoldList->firstContigInScaffold);
}

int DeleteContigList(ContigList *contigList)
{
//free memory in contigList
	free(contigList->contig);
	if (contigList->numInName <0) free(contigList->nameTable);
	else free(contigList->numToIndex);
	if (contigList->FastQDataRetained) {
		free(contigList->FastQData); //when using LoadFastQ
		contigList->FastQDataRetained=0;
	}
}



