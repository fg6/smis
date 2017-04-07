//Joe Henson
//This file carries out various operations on the scaffoldList graph, e.g. checking for juctions.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "gapest.h"
#include "readlist.h"
#include "scaffoldlist.h"
#include "refine.h"

//////////////////////prelude: some display things most of which probably shouldn't be here any more

boolmy verbose;

void DisplayContigsInScaffold(Scaffold *scaf)
{
	ContigInScaffold *cis;
	int32my i;
	for (cis = scaf->frontContig,i=0; cis !=NULL && i<100; cis=cis->backwards,++i) printf(" %s ",cis->contig->name);
}

void ValidateReadPositions(Scaffold *scaf1,Scaffold *scaf2, int dir1, int dir, int gapSize)
{
	Read *read;
	HalfEdge *he;

	for (he= dir ? (dir1?scaf1:scaf2)->frontHalfEdge : (dir1?scaf1:scaf2)->backHalfEdge; he !=NULL ; he = he->next)
	{
		for (read=he->read; read!=NULL; read=read->next)
		{
			if (read->scaffoldEnd > scaf1->len + scaf2->len + gapSize || read->scaffoldStart < 0) {
				printf("(WARNING) OffsetReadPositionInHalfEdgeList: out of bounds, dir1 %d new len %ld gap %ld start %ld end %ld len1 %ld len2 %ld \n",
					dir1, scaf1->len + scaf2->len + gapSize, gapSize, read->scaffoldStart,read->scaffoldEnd,scaf1->len,scaf2->len);
			}
		}
	}
}

void ValidateHalfEdgeList(Scaffold *scaf1, int listDir, Scaffold *scaf2, int thorough)
{
	HalfEdge *he1,*he2;
	HalfEdge **previous;
	HalfEdge *newClosest;
	HalfEdge *list1,*list2;
	int32my num1, i;
	if (scaf1==scaf2) printf("MergeHalfEdgeLists: scaffolds are identical \n \n \n ");
	if (listDir)
	{
		list1= scaf1->frontHalfEdge; //NB:heads of lists may change when halfedges are deleted.
		num1=  scaf1->frontEdgeNum;
	}
	else
	{
		list1= scaf1->backHalfEdge;
		num1=  scaf1->backEdgeNum;	
	}
	for (he1=list1,i=0; he1 !=NULL; he1=he1->next,i++)
	{
		if (he1->next != NULL && he1->next->prev != he1)
			printf("\n(WARNING) incorrect prev pointer \n \n \n ");
		if (he1->oppositeScaffold == scaf1)
			printf("\n(WARNING) halfedge points to scaf1 \n \n \n ");
		if (he1->otherEnd->oppositeScaffold != scaf1)
			printf("\n(WARNING) Opposite scaf of other end is not scaf1 \n \n \n ");
		if (thorough)
		{
			list2 = he1->otherEnd->dir ? he1->oppositeScaffold->frontHalfEdge : he1->oppositeScaffold->backHalfEdge;
			for (he2=list2; he2!=he1->otherEnd && he2 !=NULL; he2=he2->next) {}
			if (he2==NULL) printf("\n(WARNING) other end not in list of opposite Scaffold \n \n \n ");
		}
	}
	if (i != num1 ) printf("\n(WARNING) actual number of halfedge %d whereas ",i);
	if (verbose) printf("recorded number %d \n ",num1);
}

void ValidateAllHalfEdges(ScaffoldList *scaffoldList)
{
	Scaffold *scaf;

	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
	{
		ValidateHalfEdgeList(scaf, 0, NULL, 1);
		ValidateHalfEdgeList(scaf, 1, NULL, 1);
	}
}

void ReportOnEdges(ScaffoldList* scaffoldList)
{
// each edge has properties like edgeQual (==0 if bad edge), an estimated gap size etc.
// this functions resets them for all edges. Used when parameters that affect edge quality are changed.
	Scaffold *scaffold;
	HalfEdge *he;
	int32my total=0, good=0; 
	for (scaffold=scaffoldList->scaffold;scaffold!=NULL;scaffold=scaffold->next) {
		for (he=scaffold->frontHalfEdge; he != NULL ; he=he->next) {
			++ total;
			if (he->edge->numGoodPair > 5) ++ good;
		}
		for (he=scaffold->backHalfEdge; he != NULL ; he=he->next) {
			++ total;
			if (he->edge->numGoodPair > 5) ++ good;
		}
	}
	printf("ReportOnEdges: total number of halfedges: %" PRI32" known good edges: %"PRI32".\n", total/2, good/2);
}

//////////////////////////////////the useful subroutines //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

////SECTION 1///////////////Resolving cases in which >1 edge (different directions) connect two scaffolds ("multiple edges")////////////////////////////////////

int32my FindWeakMultipleEdges(HalfEdge *multi[4], float dirRatio)
{
//passed a list of 4 half edges (some of which may be null)
//between the same pairs of contigs, attempts to find the real one,
//setting the others' qualities to 0.
//returns ONE LESS than # good edges left afterwards (i.e. number of remaining conflicts)
	int32my iNum,bestNum,i,numLeft=0;
	
	bestNum=0;
	for (i=0;i<4;++i) {
		if (multi[i]==NULL) continue;
		iNum=  multi[i]->edge->numGoodPair;
		if (iNum > bestNum) bestNum=iNum;
	}
	for (i=0;i<4;++i) {
		if (multi[i]==NULL) continue;
		iNum=multi[i]->edge->numGoodPair;
		if ((float) dirRatio * iNum < bestNum) {if (verbose) printf("-R)"); multi[i]->edge->edgeQual=0; multi[i]=NULL; }
	}
	for (i=0;i<4;++i)
		if (multi[i]!=NULL) ++numLeft;
	//CAUTIOUS APROACH
	if (numLeft > 1)	{
		for (i=0;i<4 && multi[i]==NULL ;++i) ;
		if (verbose) printf("Unresolved multi edge between %s and %s\n",
			multi[i]->otherEnd->oppositeScaffold->frontContig->contig->name,
			multi[i]->oppositeScaffold->frontContig->contig->name);
	}
	if (numLeft > 1 && numLeft < 4) {
		if (verbose) printf("-U)");
		for (i=0;i<4;++i) if (multi[i]!=NULL) multi[i]->edge->isMultiEdge=1;
	}
	if (numLeft == 4) {
		if (verbose) printf("-T)");
		for (i=0;i<4;++i) if (multi[i]!=NULL) multi[i]->edge->isMultiEdge=2;
	}
	return(numLeft - 1);
}

void FindContainedScaffolds(ScaffoldList *scaffoldList, Scaffold *scaf, MergeParams *params)
{
// looks at connections in a scaffold in an attempt to find smaller scaffolds which could
// be inserted into gaps in larger ones.  If so does the insertion.

	HalfEdge *he, *he2;
	Scaffold *opScaf;
	ContigInScaffold *cis;
	int32my frontDir;
	int32my gapPos, maxGap=0;
	int32my containPairThreshold;
	int32my containTotalThreshold;
	float pos1,pos2;
	float discrepancy;

	containPairThreshold=params->edgeParams.numPairThreshold > 10 ? 10 : params->edgeParams.numPairThreshold; //more or less arbitrary choice here.
	containTotalThreshold= 3 * params->edgeParams.numPairThreshold / 2 ;
	
	if (scaf->contigNum == 1) return; //no gaps for other scafs to be contained in.
	for (cis=scaf->frontContig; cis->backwards!=NULL; cis=cis->backwards)
		if (cis->backwardGap > maxGap) maxGap = cis->backwardGap;
	if (maxGap < 220) return; //nothing will fit into gaps this small.
	//printf("FindContainedScaffolds\n");
	for (he=scaf->backHalfEdge;he!=NULL;he=he->next)
	{
		if (!he->edge->edgeQual && he->edge->numPair < containPairThreshold) continue;
		//if (!he->edge->edgeQual && he->edge->gapSize==-1) QuickEstimateGapNoCutoff(he); //not all these have an estimate already but if they have don't waste time on a more accurate one yet
		if (he->oppositeScaffold->len > maxGap + 20) continue;
		opScaf=he->oppositeScaffold;
		frontDir=1 ^ he->otherEnd->dir; // dir of an edge from opScaf to the front of scaf, if opScaf is contained
		for (he2=scaf->frontHalfEdge; he2!=NULL && !(he2->oppositeScaffold==opScaf && 
				he2->otherEnd->dir==frontDir);he2=he2->next) {}
		if (he2!=NULL && he2->edge->numPair >= containPairThreshold && 
				he->edge->numPair + he2->edge->numPair >= containTotalThreshold ) {
			//if (!he2->edge->edgeQual && he2->edge->gapSize==-1) QuickEstimateGapNoCutoff(he2);
			pos1=EstimateRelativePosition(scaf, he);
			pos2=EstimateRelativePosition(scaf, he2);
			discrepancy=pos2 - pos1 - opScaf->len;
				//printf("FindContainedScaffolds: scaf %s (%s) len %d contained in %s (%s) len %d \n", 
				//	opScaf->frontContig->contig->name, opScaf->backContig->contig->name, opScaf->len,
				//	scaf->frontContig->contig->name, scaf->backContig->contig->name, scaf->len);
				//printf("he  0->%d gap %f num %d good %d bio %d range %" PRI32 " rel pos %f\n",
				//he->otherEnd->dir, he->edge->gapSize, he->edge->numPair, he->edge->numGoodPair, he->edge->numBiotin, GetHalfEdgeRangeAllReads(he), pos1);
				//printf("he2 1->%d gap %f num %d good %d bio %d range %" PRI32 " rel pos %f\ndiscrepancy %f \n",
				//he2->otherEnd->dir, he2->edge->gapSize, he2->edge->numPair, he2->edge->numGoodPair, he2->edge->numBiotin, GetHalfEdgeRangeAllReads(he2), pos2, discrepancy);
			if (abs(discrepancy) < 3000) {
				// if (verbose) printf("XXXXXXXX Small scaf %s len %"PRI32" contained in scaf %s len %"PRI32"\n",opScaf->frontContig->contig->name, opScaf->len, scaf->frontContig->contig->name, scaf->len);
				//if (verbose)printf("------- pos 1 %f pos2 %f\n",pos1,pos2);
			}
			if (abs(discrepancy) < 3000 && IsRoomInGap(scaf, &cis, &gapPos, pos1 -3000, pos2 +3000, opScaf->len + 20)) {
				if (verbose) printf("Small scaf len %"PRI32" fits in gap in scaf len %"PRI32" at position %"PRI32"\n",opScaf->len,scaf->len, gapPos);
				if (verbose) printf("------- pos 1 %f pos2 %f\n",pos1,pos2);
				if (he->edge->isMultiEdge==2) 
					{if (verbose) printf("blocked by direction uncertainty (multiedge) \n");}
				else InsertScaffold(scaffoldList, scaf, opScaf, he, he2, cis, gapPos, params);
			}
		}
	}	
}

int32my ResolveMultipleEdgesFromScaffold(ScaffoldList *scaffoldList, Scaffold *scaf, MergeParams *params)
{
// often there is >1 edge between the same scaffolds
// with different directional properties.  Only one can be correct and we don't want the spurious ones
// to get interfere with anything.  If we are sure which one is the correct one we want to set the others edgeQual=0,
// otherwise we mark all remaining ones so that they can block merges but not be merged themselves.
// this routine does this for edges from one scaffold (these need to be redone after a merge).
	int32my bestNum,iNum,bestEdge,i,j;
	HalfEdge *multi[4], *he, *he2;
	Scaffold *opScaf;
	int32my numUnresolved=0;

	//first reset flags on all edges which are about to be checked and marked/resolved.
	for (he=scaf->backHalfEdge;he!=NULL;he=he->next)  he->edge->isMultiEdge=0;
	for (he=scaf->frontHalfEdge;he!=NULL;he=he->next) he->edge->isMultiEdge=0;
	//then find cases of multiple edges and pass them on to be resolved or flagged.
	for (he=scaf->backHalfEdge;he!=NULL;he=he->next)
	{
		if (!he->edge->edgeQual) continue;
		for (i=0;i<4;++i) multi[i]=NULL;
		opScaf=he->oppositeScaffold;
		for (he2=he->next; he2!=NULL && he2->oppositeScaffold!=opScaf;he2=he2->next) {}
		if (he2!=NULL && he2->edge->edgeQual)
			if (!he->otherEnd->dir) { multi[0]=he; multi[1]=he2;} else { multi[1]=he; multi[0]=he2;}
		for (he2=scaf->frontHalfEdge; he2!=NULL && he2->oppositeScaffold!=opScaf;he2=he2->next) {}
		if (he2!=NULL && he2->edge->edgeQual) {
			if (!he->otherEnd->dir) multi[0]=he; else multi[1]=he;
			if (!he2->otherEnd->dir) multi[2]=he2; else multi[3]=he2;
			for (he2=he2->next; he2!=NULL && he2->oppositeScaffold!=opScaf;he2=he2->next) {}
			if (he2!=NULL)
				if (!he2->otherEnd->dir) multi[2]=he2; else multi[3]=he2;
		}
		/*if (scaf->len > opScaf->len) {
			if (multi[1]!=NULL && multi[2]!=NULL) {
				if (multi[0]==NULL !! multi[4]==NULL)
					inserted=TryInsert(scaf,opScaf, multi[1], multi[2]);
			} else {
                        	if (multi[0]!=NULL && multi[4]!=NULL)
					inserted=TryInsert(scaf,opScaf, multi[0], multi[4]);
			}
		}*/
		if ( /* !inserted && */(multi[0]!=NULL || multi[1]!=NULL)) {if (verbose) printf("(M-");numUnresolved+=FindWeakMultipleEdges(multi, params->dirRatio);}
	}
	//we have covered all muliple edges involving back edges but not those involving only front edges.
	for (he=scaf->frontHalfEdge;he!=NULL;he=he->next)
	{
		if (!he->edge->edgeQual) continue;
		for (i=0;i<4;++i) multi[i]=NULL;
		opScaf=he->oppositeScaffold;
		for (he2=he->next;he2!=NULL && he2->oppositeScaffold!=opScaf;he2=he2->next) {}
		if (he2!=NULL && he2->edge->edgeQual) { multi[0]=he; multi[3]=he2;}
		if (multi[0]!=NULL) {if (verbose) printf("(M)");numUnresolved+=FindWeakMultipleEdges(multi, params->dirRatio);}
	}
	//printf("\n");
	return(numUnresolved > 0);
}


void ResolveMultipleEdges(ScaffoldList *scaffoldList, MergeParams *params)
{
// often there is >1 edge between the same scaffolds
// with different directional properties.  Only one can be correct and we don't want the spurious ones
// to get interfere with anything.  If we are sure which one is the correct one we want to set the others edgeQual=0,
// otherwise we mark all remaining ones so that they can block merges but not be merged themselves.
// this routine does this for edges from all scaffolds (these need to be redone at the beginning of each stage I think).
	Scaffold *scaf;
	int32my numMulti=0;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next) {
		numMulti+=ResolveMultipleEdgesFromScaffold(scaffoldList, scaf, params);
		FindContainedScaffolds(scaffoldList,scaf, params);
	}
	if (verbose) printf("(NUM-MULTI) there are %ld unresolved multiple edges \n", numMulti);
}

void FindAllContainedScaffolds(ScaffoldList *scaffoldList, MergeParams *params)
{
// often there is >1 edge between the same scaffolds
// sometimes this is because one scaffold belong inside another one.
// this routine detects this and does the insertions.
	Scaffold *scaf;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=scaf->next)
		FindContainedScaffolds(scaffoldList,scaf, params);
}

////SECTION 2///////////////merging scaffolds//////////////////////////////////////////////////////

void FlipHalfEdgeList(Scaffold *scaf, int dir)
{
//flip direction of edges and correct read's position wrt scaffold now it is flipped
	HalfEdge *he;
	Read *read;
	int32my swap;
	for (he= dir ? scaf->frontHalfEdge : scaf->backHalfEdge; he !=NULL ; he = he->next)  {
		he->dir ^= 1; 
		for (read=he->read; read!=NULL; read=read->next)  {
			read->dir_s^= 1;
			swap=read->scaffoldStart;
			read->scaffoldStart=scaf->len - read->scaffoldEnd;
			read->scaffoldEnd=scaf->len - swap;
		}
	}
}

void FlipScaffold(Scaffold *scaf)
{
// flips all the directional information in scaffold
// (direction of ContigInScaffold lists, halfedgelists, read position and direction on the scaffold)
	HalfEdge *he,*heSwap;
	ContigInScaffold *cSwap,*con;
	int32my swap;
	char cswap;
	//first reverse end points
	cSwap=scaf->frontContig;
	scaf->frontContig=scaf->backContig;
	scaf->backContig=cSwap;	
	//these lists go forwards to front and backwards to back so start on what was the front and go to what used to be back
	for (con=cSwap; con!=NULL; con=con->forwards)
	{
		con->dir ^= 1;
		cSwap= con->forwards;
		con->forwards= con->backwards;
		con->backwards=cSwap;
		swap=con->forwardGap;
		con->forwardGap=con->backwardGap;
		con->backwardGap=swap;
		cswap=con->forwardMark;
		con->forwardMark=con->backwardMark;
		con->backwardMark=cswap;
		heSwap=con->forwardHe;
		con->forwardHe=con->backwardHe;
		con->backwardHe=heSwap;
	}
	// swap front and back halfEdge lists
	heSwap=scaf->frontHalfEdge;
	scaf->frontHalfEdge=scaf->backHalfEdge;
	scaf->backHalfEdge=heSwap;
	// swap the count of front and back edges
	swap=scaf->frontEdgeNum;
	scaf->frontEdgeNum=scaf->backEdgeNum;
	scaf->backEdgeNum=swap;
	//flip direction of edges and correct read's position wrt scaffold now it is flipped
	FlipHalfEdgeList(scaf, 0);
	FlipHalfEdgeList(scaf, 1);
}

boolmy_fast ShouldMerge ( HalfEdge *he1, HalfEdge *he2 )
// true if the two half edges point to same Scaffold and have same
// directional properties
{
	return (he1->oppositeScaffold == he2->oppositeScaffold && 
		he1->dir == he2->dir && 
		he1->otherEnd->dir == he2->otherEnd->dir);
}

void MergeEdges(HalfEdge *he1,HalfEdge *he2)
// merge the edges given by the half edges into that associated to he1 
// (does nto remove he2 from its list)
{
	HalfEdge *other1, *other2; //the twin hafledgees to he1 and he2
	Read *read,*otherRead;

	other1=he1->otherEnd;
	other2=he2->otherEnd;
	if (he1==he2) {
		printf("\n(WARNING) MergeEdges: halfedges are identical \n \n \n "); 
		exit(0);
	}
	//in the halfedges themselves only the list of reads need to be combined
	he1->numReads += he2->numReads; //SetEdge needs this to be updated
	other1->numReads += other2->numReads;
	//NB: all readlists have at least one member
	for (read=he2->read, otherRead=other2->read; read->next!=NULL; read=read->next, otherRead=otherRead->next) {} //get ends of lists
	read->next=he1->read; //stick he1's list to end
	he1->read=he2->read; //put this combined list as he1's list
	otherRead->next=other1->read;
	other1->read=other2->read;
}

void MergeHalfEdgeLists (Scaffold *scaf1, int listDir, Scaffold *scaf2)
{
//merge the front or back halfedge lists into scaf1's list (deletes edges between the two scaffolds also).
//listDir= are we erging front (1) or back (0) halfedge lists of the scafs?

// (NB: two merging strategies were (a) to calc all closests from all scafs first, and then look for ones that agreed in both directions
// and merge, or (b) calcing only the closest and its reverse, merging and then recalcing the closests
// this way is not only better at finding inconsistent edges and way faster but also means that closest pointers
// do not have to be updated on merge, which led to all kinds of hideous complications here.	
	HalfEdge *he1,*he2;
	HalfEdge **list1ptr,**list2ptr;
	int *num1Ptr, *num2Ptr;

	if (listDir) {
		list1ptr= &(scaf1->frontHalfEdge); //NB:heads of lists may change when halfedges are deleted.
		list2ptr= &(scaf2->frontHalfEdge);
		num1Ptr=  &(scaf1->frontEdgeNum);
		num2Ptr=  &(scaf2->frontEdgeNum);
	} else {
		list1ptr= &(scaf1->backHalfEdge);
		list2ptr= &(scaf2->backHalfEdge);
		num1Ptr=  &(scaf1->backEdgeNum);
		num2Ptr=  &(scaf2->backEdgeNum);		
	}

	if (*list2ptr == NULL) {//if list 2 is empty of half edges, no merge, only delete edges between the two merging scaffolds
		if (*list1ptr==NULL)
			for (he1=*list1ptr; he1 !=NULL; he1=he1->next) { 
				if (he1->oppositeScaffold==scaf2) RemoveEdge(he1); //remember any opposite scaffold==scaf2 has been changed to scaf1
				if (he1->oppositeScaffold==scaf1) printf("(WARNING) MergeHalfEdgeLists: reflexive edge \n \n \n ");
			}
		return;
	}
	// otherwise, edges in the lists to the same scaf and with the same directions need to be merged, other edges from scaf2 need to be added to scaf1's list
	// edges between the two contigs must be deleted, 
	// and all pointers to scaf2 in half edges (which are all of form he->otherEnd->oppositeScaffold where he is from scaf2) need to be replaced with scaf1
	for (he2=*list2ptr; he2 != NULL ; he2=he2->next) //first get rid of the edges from scaf2 to scaf1
		if (he2->oppositeScaffold==scaf1) RemoveEdge(he2);

	for (he1=*list1ptr; he1 !=NULL; he1=he1->next)	{
		//get rid of the edges from scaf1 to scaf2
		if (he1->oppositeScaffold==scaf2) { 
			RemoveEdge(he1);
			continue;
		}
		//and merge the others if there is a matching edge from scaf2.
		for (he2=*list2ptr; he2 !=NULL && ! ShouldMerge(he1,he2); he2=he2->next) ;
		if (he2!=NULL) { //this means that a match to he1 from scaf1 was found in list from scaf 2
			MergeEdges(he1,he2);
			RemoveEdge(he2);
		}
	}
	// edge that remain in (i.e. have not been merged and deleted from) list 2 will be added to scaf1's list 1
	// but their otherEnd's have the wrong oppositescaffold, it should now be scaf1 not scaf2.  Fix this first.
	for (he2=*list2ptr; he2 != NULL ; he2=he2->next) he2->otherEnd->oppositeScaffold=scaf1;
	*num1Ptr+=*num2Ptr;
	if (*list1ptr == NULL)
		*list1ptr=*list2ptr;
	else {
		for (he1=*list1ptr; he1->next !=NULL ; he1=he1->next) {}
		he1->next=*list2ptr;
		if (*list2ptr !=NULL) (*list2ptr)->prev=he1;
	}

	ValidateHalfEdgeList(scaf1, listDir, scaf2, 0);
}

void OffsetReadPositionInHalfEdgeList(Scaffold *scaf, int dir, long offset)
{
//correct reads' positions wrt scaffold
	HalfEdge *he;
	Read *read;
	for (he= dir ? scaf->frontHalfEdge : scaf->backHalfEdge; he !=NULL ; he = he->next)
		for (read=he->read; read!=NULL; read=read->next) {
			read->scaffoldStart+=offset;
			read->scaffoldEnd+=offset;
		}
}

void UpdateReads(Scaffold *scaf1,Scaffold *scaf2, int dir1, int gapSize)
// updates the position of reads on scaffolds
// in the case that scaffold1 and 2 are being merged
{
	if (dir1) {
		OffsetReadPositionInHalfEdgeList(scaf1, 0, scaf2->len + gapSize);
		OffsetReadPositionInHalfEdgeList(scaf1, 1, scaf2->len + gapSize);
	} else {
		OffsetReadPositionInHalfEdgeList(scaf2, 0, scaf1->len + gapSize);
		OffsetReadPositionInHalfEdgeList(scaf2, 1, scaf1->len + gapSize);
	}
}

void UpdateReadsForInsertion(Scaffold *little, int32my pos)
// updates the position of reads on scaffolds
// in the case that slittle is inserted into big at position pos
{
	OffsetReadPositionInHalfEdgeList(little, 0, pos);
	OffsetReadPositionInHalfEdgeList(little, 1, pos);
}

void ResetEdgesFromScaffold(Scaffold *scaf, SetEdgeParams *params)
{
// recalculates gap size and quality for all edges from scaf.
	HalfEdge *he;
	for (he=scaf->frontHalfEdge; he != NULL ; he=he->next) SetEdge(he,he->edge,params);
	for (he=scaf->backHalfEdge ; he != NULL ; he=he->next) SetEdge(he,he->edge,params);
	//PutBadHalfEdgesAtBackofLists(scaf); //reorders edges so that edgeQual==0 edges are at the back. Also breaks everything.
}

int32my GapMLE(Scaffold *scaf, HalfEdge *he)
{
	int32my *overhangs, *thisOverhang, len1, len2;
	Distribution **whichDist, **thisDist, **distList; // whichDist gives pointer to appropriate dist for each sample
	int32my newGap,i;
	int32my numLibs;
	Read *read, *read2;
	PairLib *firstLib;
	//int32my newGap;
	Window window;

	firstLib=he->read->lib - he->read->lib->libIndex;  //NTS: the data structures here are a tangled ball of string
	//  just can't be bothered to change it or refer to the readlist at this point here when it's not necessary anywhere else.
	numLibs=firstLib->numLibs;
	len1= scaf->len;
	len2= he->oppositeScaffold->len;
	window=GetWindow(len1, len2);

	overhangs= MallocWithAssert(sizeof(int32my) * he->numReads,"GapMLE, overhangs");
	whichDist= MallocWithAssert(sizeof(Distribution *) * he->numReads,"GapMLE, whichDist");
	distList=  MallocWithAssert(sizeof(Distribution *) * numLibs,"GapMLE, distList");
	
	for (i=0, thisDist=distList;  i < numLibs;  ++i, ++thisDist) //creating list of dists here
		*thisDist=&(firstLib[i].dist); //that should do it.

	for(read=he->read, thisOverhang=overhangs, thisDist= whichDist; read!=NULL;
			 read=read->next, ++thisOverhang, ++ thisDist) {
		read2=GetMate(read);
		*thisOverhang = (read->dir_s ?  read->scaffoldStart : (len1 - read->scaffoldStart))
				+ (read2->dir_s ?  read2->scaffoldStart : (len2 - read2->scaffoldStart)) ;
		*thisDist = &(read->lib->dist);
		//printf("g(%"PRI32")=%"PRI32" [%"PRI32",%"PRI32"]", read->lib->libIndex, read->lib->size - *thisOverhang,read->scaffoldStart,read2->scaffoldStart);
	}

	if (he->edge->gapSize==-1) QuickEstimateGapNoCutoff(he);
		
	newGap = MultiDistributionMLE(overhangs, he->numReads, distList, numLibs, whichDist, window, he->edge->gapSize-2000.0, he->edge->gapSize+2000.0, 400);
	newGap = MultiDistributionMLE(overhangs, he->numReads, distList, numLibs, whichDist, window, newGap-400.0, newGap+400.0, 80);
	newGap = MultiDistributionMLE(overhangs, he->numReads, distList, numLibs, whichDist, window, newGap-80.0, newGap+80.0, 20); 
	//newGap = MultiDistributionMLE(overhangs, he->numReads, distList, numLibs, whichDist, window, newGap-10.0, newGap+10.0, 1);

	/*newGap = MaxLikelihood(&(lib->dist), overhangs, he->numReads, window, he->edge->gapSize-1000.0, 
		he->edge->gapSize+1000.0, 200);
	newGap = MaxLikelihood(&(lib->dist), overhangs, he->numReads, window, newGap-200.0, 
		newGap+200.0, 20);
	newGap = MaxLikelihood(&(lib->dist), overhangs, he->numReads, window, newGap-20.0, 
		newGap+20.0, 1);*/
	if (verbose) printf("GapMLE: readNum=%d, previous gap= %f new gap= %d\n",he->numReads, he->edge->gapSize, newGap);
	//if (newGap > 0 && he->edge->gapSize < newGap) printf("GapMLE: jjjjjjjjjj newGap larger than old\n");
	//if (newGap > 0 && he->edge->gapSize < newGap - 1000) printf("GapMLE: ssssssssss newGap much larger than old\n");
	free(overhangs);
	free(whichDist);
	free(distList);
	return (newGap);
}

int InsertScaffold(ScaffoldList* scaffoldList, Scaffold *bigScaf, Scaffold *littleScaf, HalfEdge *he1, HalfEdge *he2, ContigInScaffold *backwardInsertCis, int32my gapStartPos, MergeParams *params)
// inserts littleScaf into bigScaf before insertCis which has position gapPos relative to the start of bigScaf.
// he1 is the back-directed edge from big to small and he2 is the front one.
// at the moment the littleScaf must be at least 20bp longer than the gap into which it is put.
//returns 1 if yes equals no.
{
	ContigInScaffold *forwardInsertCis;
	int32my gap, meanPos, backwardInsertGap, forwardInsertGap;
	int32my posToInsert, gapEndPos;
	int32my estimatedStartPos, estimatedEndPos;

	if (bigScaf->contigNum < 2) {
		printf("InsertScaffold: %s contains no gaps to insert into\n", bigScaf->frontContig->contig->name);
		exit(1);
	}
	if (backwardInsertCis == NULL) {
		printf("InsertScaffold: contig in scaffold NULL\n");
		exit(1);
	}
	////////////////////////////////
	// estimate the two new gap sizes
	// make better estimate of gap size and use it to find end end of littleScaf relative to front of big
	if (!he1->edge->edgeQual && he1->edge->gapSize==-1) QuickEstimateGapNoCutoff(he1);
	if (!he2->edge->edgeQual && he1->edge->gapSize==-1) QuickEstimateGapNoCutoff(he2);
	he1->edge->gapSize=GapMLE(bigScaf, he1);
	he2->edge->gapSize=GapMLE(bigScaf, he2);
	estimatedStartPos=he1->edge->gapSize + bigScaf->len;
	estimatedEndPos= - he2->edge->gapSize;
	gap=backwardInsertCis->backwardGap; //the gap into which littleScaf is to be placed
	gapEndPos=gapStartPos+gap;

	//but we know the length of the contig for sure while these are just estimates.
	if (verbose) printf("InsertScaffold: read pairs give %"PRI32" ++ %"PRI32" for scaf len %"PRI32" with gap starting at %"PRI32" size %"PRI32"\n", estimatedStartPos, estimatedEndPos, littleScaf->len, gapStartPos, gap);
	estimatedStartPos= (estimatedStartPos + (estimatedEndPos - littleScaf->len)) / 2 ;
	estimatedEndPos= estimatedStartPos + littleScaf->len;
	if (verbose) printf("InsertScaffold: averaged %"PRI32" ++ %"PRI32" for scaf len %"PRI32" \n",estimatedStartPos,estimatedEndPos,littleScaf->len);
	//estimate where the little scaff should be

	//assuming here that the gap need not be widened, which would be a pain
	backwardInsertGap= estimatedStartPos - gapStartPos; //backwards from the back of the cis.
	forwardInsertGap= gapEndPos - estimatedEndPos; //forwards from the end of gap
	if (backwardInsertGap < 10) {
		forwardInsertGap += backwardInsertGap - 10;
		backwardInsertGap = 10;
	}
	if (forwardInsertGap  < 10) {
		backwardInsertGap += forwardInsertGap - 10;
		forwardInsertGap = 10;
	}
	if (verbose) printf("InsertScaffold: gaps +%"PRI32"+scaf+%"PRI32"+ \n",backwardInsertGap,forwardInsertGap);
	////////////////////////////////
	// now work on the scaffold to do the insertion using these values.
	forwardInsertCis= backwardInsertCis->backwards;
	//if the front directed edge from big goes to the front of little then we need to flip little.
	if (he2->otherEnd->dir) FlipScaffold(littleScaf);
	// output to track resulting join in results.
	if (verbose) printf("Inserting %s (orientation to be= %d) in scaf starting with %s\n", littleScaf->frontContig->contig->name, littleScaf->frontContig->dir, bigScaf->frontContig->contig->name);
	if (verbose) printf("To be inserted between contigs %s (%d) and %s (%d)\n",
		backwardInsertCis->contig->name,backwardInsertCis->dir, forwardInsertCis->contig->name,forwardInsertCis->dir);
	//lengths and read positions need updating
	UpdateReadsForInsertion(littleScaf, gapStartPos);
	bigScaf->contigNum += littleScaf->contigNum;
	bigScaf->numReads  += littleScaf->numReads;
	//now combine the edge lists -- this is where most of the hard work gets done.
	MergeHalfEdgeLists (bigScaf, 0, littleScaf);
	MergeHalfEdgeLists (bigScaf, 1, littleScaf);
	// it remains to update the contig in scaffold list for the combined scaffold.
	littleScaf->backContig->backwards=forwardInsertCis;
	forwardInsertCis->forwards=littleScaf->backContig;
	littleScaf->backContig->backwardHe=he2->otherEnd;
	forwardInsertCis->forwardHe=he2;
	forwardInsertCis->forwardGap=littleScaf->backContig->backwardGap=forwardInsertGap;
	forwardInsertCis->forwardMark=littleScaf->backContig->backwardMark=params->mergeMark;

	littleScaf->frontContig->forwards=backwardInsertCis;
	backwardInsertCis->backwards=littleScaf->frontContig;
	littleScaf->frontContig->forwardHe=he1->otherEnd;
	backwardInsertCis->backwardHe=he1;
	backwardInsertCis->backwardGap=littleScaf->frontContig->forwardGap=backwardInsertGap;
	backwardInsertCis->backwardMark=littleScaf->frontContig->forwardMark=params->mergeMark;

	/*if(backwardInsertGap > 40000 || forwardInsertGap > 40000) {
		printf("NO!!!\n");
		exit(1);
	}*/

	if (		backwardInsertCis->backwards == backwardInsertCis || 
			backwardInsertCis->forwards  == backwardInsertCis ||
			forwardInsertCis->forwards  == forwardInsertCis ||
			forwardInsertCis->forwards  == forwardInsertCis ||
			littleScaf->frontContig->forwards == littleScaf->frontContig ||
			littleScaf->frontContig->backwards == littleScaf->frontContig ||
			littleScaf->backContig->forwards == littleScaf->backContig ||
			littleScaf->backContig->backwards == littleScaf->backContig ) {
		printf("InsertScaffold: cis points to self\n");
		exit(1);
	}
		//it remains to take scaffold2 out of the scaffold list, and to update the edge info (quality, inconistent read pairs, gap size) from scaf1.
	RemoveScaffold(scaffoldList, littleScaf);
	ResetEdgesFromScaffold(bigScaf,& params->edgeParams);
	ResolveMultipleEdgesFromScaffold(scaffoldList, bigScaf, params);
	//edges to this scaffold (but not others apart from those) may have become multiple because of merge.  This updates the marks
	if (verbose) printf("Done insert \n");
	return(0);
}


int MergeScaffolds(ScaffoldList* scaffoldList, Scaffold *scaffold1, HalfEdge *he1, MergeParams *params)
{
//Merges the two scaffolds connected by he1 combining their information (including edges) together
//the gap between the scaffolds in the contig in scaffold list is marked with value "mark" which is saved in contig list output later.
//returns 1 if merge was prevented.
	HalfEdge *he2;
	Scaffold *scaffold2;
	ContigInScaffold *join1, *join2, *cis; //these two are the ends of the scaffolds to be next to each other in the combined scaffold
	boolmy_fast dir1, dir2;
	int32my gap;
	boolmy_fast flip; // if scaffold2 needs its direction reversed	

	// make better estimate of gap size
	he1->edge->gapSize=GapMLE(scaffold1, he1);

	he2=he1->otherEnd;
	dir1=he1->dir;
	dir2=he2->dir;
	gap=he1->edge->gapSize;

	/*if(gap > 40000) {
		printf("NO!!!\n");
		exit(1);
	}*/

	flip= 1 ^ dir1 ^ dir2; //if the reads have same direction their scaffolds are in different orientation, so one (scaf2) needs to be reversed
	scaffold2=he1->oppositeScaffold;
	// output to track resulting join in results.
	if (verbose) printf("-------- Merging %s with %s\n",scaffold1->frontContig->contig->name, scaffold2->frontContig->contig->name);
	if (verbose) printf("End contigs to be next to each other in new scaffold : ec %s ec %s\n",
		dir1 ? scaffold1->frontContig->contig->name: scaffold1->backContig->contig->name,
		dir2 ? scaffold2->frontContig->contig->name: scaffold2->backContig->contig->name);
	// going to merge into scaffold1 and delete 2
	// prepare to add contiginScaffolds in 2 to list in 1 -- if the reads point in the same direction need to "flip" 2
	if (flip) FlipScaffold(scaffold2);
	//lengths and read positions need updating
	//NB: here we have a choice of using the real gap size or replacing negative gap sizes with 10.
	//As negative gap sizes are not usually errors but instead a sign of a diploid region, use them
	UpdateReads(scaffold1, scaffold2, dir1, gap ) ;
	scaffold1->contigNum += scaffold2->contigNum;
	scaffold1->numReads  += scaffold2->numReads;
	//if we take negative gap sizes seriously (which we usually should) the updated length has to take this into account.
	//the alternative is scaffold1->len=scaffold1->len + scaffold2->len + ( gap > 0 ? gap : DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP );
	scaffold1->len= scaffold2->len + gap < 0 ? scaffold1->len : scaffold1->len + gap < 0 ? scaffold2->len : scaffold1->len + scaffold2->len + gap;
	//now combine the edge lists -- this is where most of the hard work gets done.
	MergeHalfEdgeLists (scaffold1, 0, scaffold2);
	MergeHalfEdgeLists (scaffold1, 1, scaffold2);
	// it remains to update the contig in scaffold list for the combined scaffold.
	if (scaffold2->len > (- gap)) { // new feature: don't add scaf2 contigs to list if it is contained in scaf1 (scaf1 should be longer than scaf2).
		if (dir1) {
			join1= scaffold1->frontContig;
			join2= scaffold2->backContig;
			if (join1 ==  join2)	{ 
				printf("\n(WARNING)  MergeScaffolds: error, contig to be joined to itself \n"); return(1);}
			join1->forwards=  join2;
			join2->backwards= join1;
			join1->forwardGap=join2->backwardGap=gap;
			join1->forwardMark=join2->backwardMark=params->mergeMark;
			join1->forwardHe=he1;join2->backwardHe=he2;
			scaffold1->frontContig=scaffold2->frontContig;
		} else {
			join1= scaffold1->backContig;
			join2= scaffold2->frontContig;
			//printf(" MergeScaffolds: join1 %s join 2 %s \n", join1->contig->name, join2->contig->name);
			if (join1 ==  join2)	{ 
				printf("\n(WARNING)  MergeScaffolds: error, contig to be joined to itself \n"); return(1);}
			join1->backwards= join2;
			join2->forwards=  join1;
			join1->backwardGap=join2->forwardGap=gap;
			join1->backwardMark=join2->forwardMark=params->mergeMark;
			join1->backwardHe=he1;join2->forwardHe=he2;
			scaffold1->backContig=scaffold2->backContig;
		}
	}

	//it remains to take scaffold2 out of the scaffold list, and to update the edge info (quality, inconistent read pairs, gap size) from scaf1.
	RemoveScaffold(scaffoldList, scaffold2);
	ResetEdgesFromScaffold(scaffold1,& params->edgeParams);
	ResolveMultipleEdgesFromScaffold(scaffoldList, scaffold1, params);
	// FindContainedScaffolds(scaffoldList,scaffold1, params); //bit of a waste of time.
	//edges to this scaffold (but not others apart from those) may have become multiple because of merge.  This updates the marks
	if (verbose) printf("Done merge in dir %d \n",dir1);
	return(0);
}


////SECTION 3///////////////finding "closests" (candidates to merge to given scaffold)/////////////////////////////////////////////////////

int FindClosest(ScaffoldList *scaffoldList, Scaffold *scaf, int dir, HalfEdge *closest, MergeParams *params)
{
//returns 0 for if closest found (and supplies halfedge in scaf->frontClosest or backClosest, 1 if no good connections, 
//2 if there is a conflict, 3 if blocked by multiedge.

	HalfEdge *he, *heList,**heClosest, *nextHe, *betweenHe;
	int32my closest_b,closest_e; //distance from relevant end of scaffold to beg/end of closest scaffold
	int32my he_b; //relative position beginning and ends of current scaffold
	int16f goodEdgeNum, conflictNum; //number of edges passing quality standard
	float closestStrength,nextStrength,thisStrength;
	float uncertainty,closestUncertainty; //amount by which reltaive positions of contigs are allowed to be shifted to account for error

	//some verbose output
	ContigInScaffold *endCis;
	endCis= dir ? scaf->frontContig : scaf->backContig;
	if (verbose) printf(">>>>>>> FindClosest: scaffold %s (len %d) in dir %d (end contig %s) \n", scaf->frontContig->contig->name, scaf->len,dir, endCis->contig->name);

	heList= dir ? scaf->frontHalfEdge : scaf->backHalfEdge;
	heClosest = dir ? &(scaf->frontClosest) : &(scaf->backClosest);
	if (dir) scaf->isFrontJunction= 0; else scaf->isBackJunction= 0; //reset this flag before seing if it is now a junction

	//find the candidate closest if it has not been provided
	//initially it is the edge with the smallest gap
	if (verbose) printf("closest candidates have gap sizes ");
	if (closest==NULL)
	{
		goodEdgeNum= 0;
		closest_b=10000000; // will the relative position of the closest scaffold to the relevant end of scaf (i.e. the gap size)
		for (he=heList; he!=NULL ; he=he->next) {
			//we will consider for joining those edges that are good quality, not to dormant scaffolds,
			//and not to scaffolds that fall entirely within the "skip zone".
			/*if (verbose) printf("he is QQ %s WQ %s len %d gap %f num %d bio %d range %d\n", 
		he->oppositeScaffold->frontContig->contig->name, he->oppositeScaffold->backContig->contig->name,
		he->oppositeScaffold->len, he->edge->gapSize,
		he->edge->numGoodPair, he->edge->numBiotin, GetHalfEdgeRange(he));*/
			if (! he->edge->edgeQual || he->oppositeScaffold->isDormant ||
				he->edge->gapSize + he->oppositeScaffold->len < params->skipLen ) continue;
			++goodEdgeNum;
			if (verbose) printf("%f,",he->edge->gapSize);
			if (he->edge->gapSize < closest_b) {
				closest=he;
				closest_b=he->edge->gapSize; //note this allows negative gap sizes to win
			}
		}
		if(! goodEdgeNum) { //if there are no good connections
			if (verbose) printf("\n<<<<<<<<1 No good edges from this scaf in this dir.\n");
			*heClosest=NULL;
			return(1);
		}
	} else {
		closest_b=closest->edge->gapSize;
	}
	if (verbose) printf("\n");
	// for the purposes of overlap checking we set the relative position of the far end of the closest scaffold
	closest_e= closest_b + closest->oppositeScaffold->len;
	if (verbose) printf("closest is %s (%s) dir %" PRIBOOLF " len %" PRI32 " gap %f num %d bio %d range %d score %f\n", 
		closest->oppositeScaffold->frontContig->contig->name, closest->oppositeScaffold->backContig->contig->name,
		closest->otherEnd->dir, closest->oppositeScaffold->len, closest->edge->gapSize,
		closest->edge->numGoodPair, closest->edge->numBiotin, GetHalfEdgeRange(closest), closest->edge->meanScore);
	//then find conflicts
	closestStrength= (float) closest->edge->numGoodPair; // GetStrengthC(scaf, heList, dir, closest);
	closestUncertainty = params->findSigmas / closest->edge->w > params->minUncertainty ? params->findSigmas / closest->edge->w : params->minUncertainty ;
	conflictNum=0;
	nextStrength=0;
	for (he=heList; he!=NULL ; he=he->next)	{
		//the same test as before on edges that should be considered for joining, and also ignore the closest
		if (!he->edge->conflictCheck || he->oppositeScaffold->isDormant || he==closest ||
			he->edge->gapSize + he->oppositeScaffold->len < params->skipLen)  continue;
		he_b = he->edge->gapSize; //the beginning of the scaffold connected by he to scaf in co-ords relative to end of scaf.
		uncertainty = params->findSigmas / he->edge->w > params->minUncertainty ? params->findSigmas / he->edge->w : params->minUncertainty ;
		// crucial lines: checks for conflicts between closest and other he's.  The definition of conflict is one of the main parts of algorithm.
		// if there is no edge (betweenHe) from closest to he, in the relevant direction, this is a conflict.
		// Also, if there is an overlap between closest and he in terms of their relative positions, this is a conflict,
		// apart from if the gap size between closest and he is negative.  Then we will allow an overlap smaller than (-ve) this gap size.
		betweenHe = GetHalfEdge(closest->oppositeScaffold, 1 ^ closest->otherEnd->dir, he->oppositeScaffold, he->otherEnd->dir);
		if ( betweenHe == NULL  ||
			he_b + uncertainty  <= closest_e + (betweenHe->edge->gapSize < 0 ? betweenHe->edge->gapSize : 0 ) - closestUncertainty )
		{
			if (verbose) printf("conflict with %s (%s) dir %" PRIBOOLF " len %" PRI32 " gap %f num %d bio %d range %d score %f\n",
				he->oppositeScaffold->frontContig->contig->name, he->oppositeScaffold->backContig->contig->name, he->otherEnd->dir,
				he->oppositeScaffold->len, he->edge->gapSize, he->edge->numGoodPair, he->edge->numBiotin, GetHalfEdgeRange(he), he->edge->meanScore);
			//the overlap may be longer than the whole scaffold -- special rule for this:
			if ( betweenHe != NULL && (- betweenHe->edge->gapSize) > he->oppositeScaffold->len ) {
				if (verbose) printf("contained in closest\n"); // here there is an edge confirming that he should be inside closest
				//NB: this will not catch all containments -- large ones will have no edge as they will be ruled out as "chimeric".
			} else {
				++ conflictNum; //+1 conflicts if passes test
				if ( (thisStrength=(float) he->edge->numGoodPair) > nextStrength) {
					nextStrength = thisStrength;
					nextHe=he;
				}
			}
		}
	}
	if (! conflictNum) { //if there is no conflict use this closest
		if (closest->edge->isMultiEdge) {if (verbose) printf("<<<<<<<< 3 closest is multiedge, cannot be used \n"); *heClosest=NULL; return(3);}
		*heClosest=closest;
		if (verbose) printf("<<<<<<<< 0 no conflict, closest found \n");
		return(0);
	}
	//otherwise do further tests on conflicts.	
	//these lines resolve remaining conflicts in roughly the SSPACE style
	if (verbose) printf("strength: %f nextstrength %f ratio %f to compare %f\n", closestStrength, nextStrength, params->strengthRatio, closestStrength * params->strengthRatio );
	if (nextStrength < closestStrength * params->strengthRatio || 
			(! nextHe->edge->edgeQual && nextStrength < closestStrength * params->edgeParams.StrengthRatioForWeakEdges))	{
		if (closest->edge->isMultiEdge) {
			if (verbose) printf("<<<<<<<< 3 closest beat strength test but is multiedge, cannot be used \n");
			*heClosest=NULL;
			 return(3);
		}
		if (verbose) printf("<<<<<<<< 0 closest beat strength test, confirmed \n");
		*heClosest=closest;
		return(0);
	}
	if ((float) nextStrength * params->strengthRatio > closestStrength ) {
		if (verbose) printf("Closest failed strength test. recalculating closest.\n <<<<<<<<");
		return (FindClosest(scaffoldList, scaf, dir, nextHe, params));
	}
	// if all these fail we have a junction
	if (dir) scaf->isFrontJunction=1; else  scaf->isBackJunction=1;	
	if (verbose) printf("<<<<<<<< 2 unresolved conflict \n");
	*heClosest=NULL;
	return(2);
}

int FindClosestAndConverseClosestNoSkip(ScaffoldList *scaffoldList, Scaffold *scaf, int dir, MergeParams *params)
{
// attempts to find the half edge to the scaffold that is closest to scaf in relative distance according to the mate pair data
// in direction specified, using parameters specified.
// if this can be done returns 1 and puts pointer to the half edge in scaf->front(back)Closest
	int rvalue;
	if (verbose) printf("FindClosestAndConverseClosestNoSkip \n");
	if (rvalue=FindClosest(scaffoldList, scaf, dir, NULL, params)) {
		//if (verbose) printf("Return value %d \n", rvalue);
		return(0); //findclosest returns 0 when a closest is found.
	}
	if (verbose) printf("Checking opposite closest for consistency { \n");
	HalfEdge *closest;
	closest= dir ? scaf->frontClosest : scaf->backClosest; //NTS: there was a point at which this was made sense... is it now now pointless to store closests?
	if (FindClosest(scaffoldList, closest->oppositeScaffold, closest->otherEnd->dir, NULL, params) || 
		(closest->otherEnd->dir ? closest->oppositeScaffold->frontClosest : closest->oppositeScaffold->backClosest) != closest->otherEnd ) {
		if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL; //don't want calling function using the bad closest pointer.
		if (verbose) printf("} opposite closest NOT consistent \n");
		return(0);
	}
	if (verbose) printf("} opposite closest consistent \n");
	if (params->preventTwoSmallScaffoldMergeThreshold > scaf->len && //a last check that removes joins between two small contigs
		params->preventTwoSmallScaffoldMergeThreshold > closest->oppositeScaffold->len) {
		if (verbose) printf("Merge prevented by preventTwoSmallScaffoldMergeThreshold \n");
		if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL;
		if (closest->otherEnd->dir) closest->oppositeScaffold->frontClosest=NULL; else closest->oppositeScaffold->backClosest=NULL;
		return(0);
	}
	return(1);
}

int FindClosestAndConverseClosest(ScaffoldList *scaffoldList, Scaffold *scaf, int dir, MergeParams *params)
{
// attempts to find the half edge to the scaffold that is closest to scaf in relative distance according to the mate pair data
// in direction specified, using parameters specified.
// if this can be done returns 1 and puts pointer to the half edge in scaf->front(back)Closest, otherwise returns 0
	if (params->skipLen == NO_SKIP) return(FindClosestAndConverseClosestNoSkip(scaffoldList, scaf, dir, params));
	MergeParams noSkipParams;
	int result1,result2;
	HalfEdge *closest;

	//first make an attempt to find closest with no skip.
	noSkipParams=*params;
	noSkipParams.skipLen=NO_SKIP; //bit of a longwinded way of doing this I suppose.
	result1=FindClosest(scaffoldList, scaf, dir, NULL, &noSkipParams);
	if (result1 == 1) return(0); //no good connections at all
	if (!result1) { //if we found a closest try for the other direction
		closest= dir? scaf->frontClosest : scaf->backClosest;
		result2=FindClosest(scaffoldList, closest->oppositeScaffold, closest->otherEnd->dir, NULL,  &noSkipParams);
		if (result2==1) { //no good connections in backwards direction
			if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL; //don't want calling function using the bad closest pointer.
			return(0);
		}
		if ((closest->otherEnd->dir ? closest->oppositeScaffold->frontClosest : closest->oppositeScaffold->backClosest) == closest->otherEnd) {
			if (params->preventTwoSmallScaffoldMergeThreshold > scaf->len && //a last check that removes joins between two small contigs
				params->preventTwoSmallScaffoldMergeThreshold > closest->oppositeScaffold->len) {
				if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL;
				if (closest->otherEnd->dir) closest->oppositeScaffold->frontClosest=NULL; else closest->oppositeScaffold->backClosest=NULL;
				return(0);
			}
			return(1);
		}
	} //that leaves cases where there was a conflict in either direction or where the closests didn't match -- couldn't find closest without skip.
	if (FindClosest(scaffoldList, scaf, dir, NULL, params)) return(0); //try with skip, true if no closest found even with skip.
	closest= dir? scaf->frontClosest : scaf->backClosest;
	if (FindClosest(scaffoldList, closest->oppositeScaffold, closest->otherEnd->dir, NULL, params) || 
		(closest->otherEnd->dir ? closest->oppositeScaffold->frontClosest : closest->oppositeScaffold->backClosest) != closest->otherEnd ) {
		if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL; //don't want calling function using the bad closest pointer.
		return(0);
	}
	if (params->preventTwoSmallScaffoldMergeThreshold > scaf->len && //a last check that removes joins between two small contigs
		params->preventTwoSmallScaffoldMergeThreshold > closest->oppositeScaffold->len) {
		if (dir) scaf->frontClosest=NULL; else  scaf->backClosest=NULL;
		if (closest->otherEnd->dir) closest->oppositeScaffold->frontClosest=NULL; else closest->oppositeScaffold->backClosest=NULL;
		return(0);
	}
	return(1);
}

////SECTION 4///////////////Overall Control of merge process///////////////////////////////////////////////////////////////


void ResetEdges(ScaffoldList* scaffoldList,SetEdgeParams *params)
{
// each edge has properties like edgeQual (==0 if bad edge), an estimated gap size etc.
// this functions resets them for all edges. Used when parameters that affect edge quality are changed.
	Scaffold *scaffold;
	HalfEdge *he;
	for (scaffold=scaffoldList->scaffold;scaffold!=NULL;scaffold=scaffold->next) {
		for (he=scaffold->frontHalfEdge; he != NULL ; he=he->next)
			if (he->oppositeScaffold > scaffold) { // the edge info is the same in both directions, but each direction has a half edge.  This selects one.
				CrossBiotinFlag(he, params->crossBiotinSWDeficit);
				LowMatchScoreFlag(he, params->matchScoreThreshold);
				SetEdge(he,he->edge,params);
			}
		for (he=scaffold->backHalfEdge; he != NULL ; he=he->next)
			if (he->oppositeScaffold > scaffold) { // the edge info is the same in both directions, but each direction has a half edge.  This selects one.
				CrossBiotinFlag(he, params->crossBiotinSWDeficit);
				LowMatchScoreFlag(he, params->matchScoreThreshold);
				SetEdge(he,he->edge,params);
			}
		//PutBadHalfEdgesAtBackofLists(scaffold);
	}
}

void FindAndMergeGoodConnections(ScaffoldList *scaffoldList, MergeParams *params)
{
// Finds the supposed closest scaffolds to some scaffold and merges the pair on some conditions
// while doing so, trys to find and merge the strongest connections (by number of confirming pairs) first.
	Scaffold *scaf, *nextScaf, *nextScafInList, *thisScaf;  //scaf is the loop index
	boolmy_fast dir, furtherDir;
	HalfEdge *closest, *oldClosest;
	Edge *frontEdge, *backEdge, *edge1, *edge2;
	int32my strength1,strength2;

	nextScafInList=scaffoldList->scaffold->next;
	for (scaf=scaffoldList->scaffold;scaf!=NULL;scaf=nextScaf) { //next scaf will be either the newly merged scaf or the next in the list if there is one
		//this first part takes the current scaf and attempts to find the strongest good "closest" connection (if there is any) at one of its two ends
		//to be merged an edge X->Y must be (a) a "closest" connection (b) closest also from Y->X (this relation is not always symmetric) and
		//(c) it must be the strongest such connection for both of the scaffolds involved, i.e. a local maximum of the strength..
		//(NB: there would be ways to avoid some of the findClosest calls here after strength tests but I doubt much time is lost).
		if (scaf->isDormant) {
			nextScaf=nextScafInList;
			if (nextScaf != NULL) nextScafInList=nextScaf->next; //update this for next time
			continue;
		}
		if (verbose) printf("Finding closests to next scaffold \n");
		FindClosestAndConverseClosest(scaffoldList, scaf, 0, params); //NB these wipe closest pointer of A if B is closest but B is not closest to A
		FindClosestAndConverseClosest(scaffoldList, scaf, 1, params);
		//if there are no closests (i.e. no connection, or a conflict, or non-symmetric) try next scaf in list
		if (scaf->frontClosest == NULL && scaf->backClosest == NULL) {
			nextScaf=nextScafInList;
			if (nextScaf != NULL) nextScafInList=nextScaf->next; //update this for next time
			continue;
		}
		//otherwise get strongest (or only) closest
		if (scaf->frontClosest == NULL) closest=scaf->backClosest; else {
			if (scaf->backClosest == NULL) closest=scaf->frontClosest; else {
				frontEdge=scaf->frontClosest->edge;  backEdge=scaf->backClosest->edge;
				closest = frontEdge->numGoodPair > backEdge->numGoodPair ?
					 scaf->frontClosest : scaf->backClosest ;
			}
		}
		//If it is good should we merge?  No, if this closest scaf has a yet stronger closest link on the other side do that one first, and so on.
		//this way we move through the linear order of the closests and find a local max in the strength, then work outwards from there.
		thisScaf=closest->oppositeScaffold; //is the strength best on the current "closest" edge or on the closest at other end of this scaffold?
		furtherDir=1 ^ closest->otherEnd->dir; //otherEnd->dir is direction from thisScaf back to scaf, we want the opposite direction to a new scaffold.
		strength1 = closest->edge->numGoodPair;
		oldClosest=closest;
		while ( 1 )
		{
			if (verbose) printf("Found candidate to merge. Finding closest at other end of candidate scaf and testing strengths...\n");
			FindClosestAndConverseClosest(scaffoldList, thisScaf, furtherDir, params);
			closest= furtherDir ? thisScaf->frontClosest : thisScaf->backClosest; //closest at the suceeding end
			if (closest == NULL) break; //if there's nothing, we can merge oldClosest
			strength2 = closest->edge->numGoodPair;
			if (strength1 >= strength2) break; //if the one we did previously is stronger we can merge it.
			//if break did not occur this closest halfedge is now the strongest good one we have found, we need to test it.
			thisScaf=closest->oppositeScaffold;
			furtherDir=1 ^ closest->otherEnd->dir;
			strength1=strength2;
			oldClosest=closest;
		}
		if (verbose) printf("\n");
		// by the time we have broken out of the previous loop, we have two scafs that can be merged
		//merge smaller into larger to best preserve the order (this is later sorted, I guess this makes it faster).
		if (oldClosest->otherEnd->oppositeScaffold->len > thisScaf->len) 
			MergeScaffolds(scaffoldList, oldClosest->otherEnd->oppositeScaffold, oldClosest, params);
		else	
			MergeScaffolds(scaffoldList, thisScaf, oldClosest->otherEnd, params);
		nextScaf=thisScaf; //go back and start again with this scaffold.
	}
}

void CollapseEdges(ScaffoldList *scaffoldList, MergeParams *params)
{
// builds scaffolds by collapsing edges for which the two scaffolds are closest to each and there are no conflicting edges
// leaving only edges to possible junctions (and some special cases involving short contigs).
	int32my i;
	Scaffold *scaffold;
	verbose=params->verbose;
	SortScaffoldListByLength(scaffoldList); //this makes sure we are attempting merges on the biggest first, although see FindAndMergeGoodConnections.
	ResetDormant(scaffoldList); //so that nothing is dormant (so that edges to them are not used in scaffolding process) at first.
	if (params->excludeSmallContigs) MakeSmallContigsDormant(scaffoldList, params->smallContigThreshold);
	ResolveMultipleEdges(scaffoldList, params); //sometimes >1 edge between a pair of scafs, with differents dirs.
	FindAndMergeGoodConnections(scaffoldList,params); //this is where the meat is.
	for (scaffold=scaffoldList->scaffold, i=0;scaffold!=NULL;scaffold=scaffold->next, ++i) ;
	if (verbose) printf("number of remaining scaffolds: %"PRI32" \n", i);
}

void DoScaffolding(ScaffoldList *scaffoldList, MergeParams *params, int numStages)
// the scaffolding process can be repeated with different parameters.
// this routine calls the scaffolding routine multiple times with different parameter sets.
{
	int32my stage, i, j, oldi,totalLen;
	Scaffold *scaf;
	for (stage=0;stage<numStages;++stage) {
		printf("DoScaffolding:Starting stage %d of scaffolding process\n", stage+1 );
		for (scaf=scaffoldList->scaffold,i=0;scaf!=NULL;scaf=scaf->next, ++i) ; //i now = number of scafs
		oldi= i + params[stage].minMergesPerStep; //just to get into first step of loop.
		// for each "stage" there is a number of "steps".  One call to CollapseEdges will not make all possible merges with the given parameters,
		//rather, it attempts merges starting from each scaffold in turn and stops when it reaches the end of list.
		//merges that happen later can remove conflicts occurring earlier so we need to go back to the beginning 
		//of list a number of times to catch everything.  The max number of steps as well as the minimum merges made in a step before stopping are parameters.
		ResetEdges(scaffoldList,&(params[stage].edgeParams));//sets gap size, quality etc.  Needs resetting here as parameters may have changed.
		ReportOnEdges(scaffoldList);

		for (j=1; j<=params[stage].maxScaffoldingSteps && i <= oldi - params[stage].minMergesPerStep; ++j) {
			printf("DoScaffolding: stage %d step %" PRI32 " \n",stage+1,j);
			if (params[stage].minMergesPerStep < 1) printf("(WARNING) DoScaffolding: minMergesPerStep= %d<1\n",params[stage].minMergesPerStep );
			CollapseEdges(scaffoldList,params + stage);
			oldi=i;
			for (scaf=scaffoldList->scaffold,i=0;scaf!=NULL;scaf=scaf->next, ++i) ;
			printf("DoScaffolding: number of remaining scaffolds after stage  %"PRI32":  %"PRI32" \n",j,i);
		}
		//FindAllContainedScaffolds(scaffoldList, params);//fills in any contigs that fit into gaps that we have missed. Doesn't need doing that often.
	}
}
