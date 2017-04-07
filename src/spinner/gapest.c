// readlist.c
// routines for finding the distribution of insert sizes
// and using these to estimate gap sizes

#define MATCH_THRESHOLD 15
#define PI 3.14159265
#define ONE_OVER_SQRT_2PI 0.39894228
//pairs of reads with a read with match score below this are not used

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "defs.h"
#include "gapest.h"
#include "readlist.h"
#include "scaffoldlist.h"

Window GetWindow(int32my len1, int32my len2)
{
//sets the parameters for the window function appropriate for two scaffolds of the specified lengths
	Window window;
	if (len1 < len2) {
		window.len1= len1;
		window.len2= len2;
	} else {
		window.len1= len2;
		window.len2= len1;
	}
	return(window);
}

double MeanValue(Distribution *d)
//gets the mean from the distribution of weight
{
	double sum=0.0;
	double *ptr;
	int32my i;
	for (i=0,ptr=d->weight; i<d->size; ++i,++ptr)
		sum+=*ptr * i;
	return ( sum / d->norm );
}

double StdDev(Distribution *d)
//gets the standard deviation from the distribution of weight
{
	double sum=0.0;
	double sumsq=0.0;
	double *ptr;
	int32my i;
	for (i=0,ptr=d->weight; i<d->size; ++i,++ptr) {
		sum+=*ptr * i; //sometimes I wish they had used something else to signify a pointer.
		sumsq+=*ptr * i * i;
	}
	sum /= d->norm; //now sum is actually mean;
	sumsq /= d->norm;
	return ( sqrt( sumsq - sum * sum ) );
}

double MaxWeight(Distribution *d)
{
	int32my i;
	double *ptr;
	double max = 0.0;
	for (i=0,ptr=d->weight; i<d->size; ++i,++ptr)
		if (*ptr > max) max=*ptr;
	return(max);
}

Distribution GetEmptyDistribution(int32my size)
//set a distribution of size size with no weight recorded
{
	int32my i;
	Distribution dist;
	dist.weight = (double *) malloc( sizeof(double) * size); // recording the frequency of each insert size from 0 up to "size" -1
	for (i=0;i<size;++i) dist.weight[i]=0.0;
	dist.size=size;
	dist.norm=0.0;
	return(dist);
}

void SetNormalDistribution(Distribution *dist, double mean, double sigma)
//set a normal distribution
{
	int32my i;
	double normFactor,var;
	var=sigma*sigma;
	normFactor= (ONE_OVER_SQRT_2PI / sigma);
	for (i=0;i<dist->size;++i)
		dist->weight[i]= normFactor * exp ( - ((float) i - mean) * ((float) i - mean) / (2.0 * var));
	dist->norm=1.0;
}

void DeleteDistribution(Distribution dist)
//free memory for dist
{
	free(dist.weight);
}

void Normalise(Distribution *dist)
//Makes a probability distribution by normalising dist to 1
{
	int32my i;
	double *weight;
	weight=dist->weight;
	for (i=0;i<dist->size;++i) weight[i] /= dist->norm;
	dist->norm=1.0;
}

void ResetNorm(Distribution *dist)
//The norm should the the sum of the weights
//This forces that to be true (corrects floating point errors etc.)
{
	int32my i;
	dist->norm=0.0;
	for (i=0;i<dist->size;++i) dist->norm += dist->weight[i];
}

void Smooth(Distribution dist, double minBlockWeight)
// intended to get rid of zeros of distributions.
// blocks of values in the distribution with aggregate weight less than minBlockWeight are given the mean weight for the block
{
	int32my i,j,k,go;
	double blockWeight;
	double min;
	min= 1/minBlockWeight;
	//go=1;
	for (i=0; i < ( dist.size -1) ;i=j-1)
	{
		if (dist.weight[i] > minBlockWeight) {j=i+2; continue;}
		blockWeight=0;
		for (j=i;j<dist.size && blockWeight < minBlockWeight ;++j)
			blockWeight+=dist.weight[j];
		if (j<dist.size) {
			//if we have found more than minBlockWeight, average minBlockWeight over all values in block and put what's left back on last value.
			min = minBlockWeight / (double) (j - i);
			dist.weight[j-1] = min + blockWeight - minBlockWeight;
		} else {
			//if not average out what we found.
			//go=0;
			min = blockWeight / (double) (j - i);
			dist.weight[j-1]= min;
		}
		//printf("min %e \n",min);
		for (k=i;k<(j-1);++k) dist.weight[k] = min;
	}
}

double WindowWeight(Window w, long x)
{
//a trapezoid shaped function depending on the window parameters, giving 0 outside window.
//the function sums to 1.
	return ( (x <= 0 ? 0.01
		: x < w.len1 ? x
		: x < w.len2 ? w.len1
		: x < (w.len1+w.len2) ? (w.len1+w.len2) - x
		: 0.01) / (double) (w.len2) );
}

double GetWeight(Distribution dist, int32my x, int32my gap, Window w, double pseudocount)
// supplies the weight from the insert size pdf dist appropriate for "overhang"x, gap gap, and
// scaffold lengths supplied by window.  *Includes the window function and pseudocounts*.
// use GetNormWithWindow to normalise these weights to get probs.
{
	return ( ( ((x+gap)>=0 && (x+gap)<dist.size) ? dist.weight[x + gap] + pseudocount : pseudocount ) //Prob for value x + gap for the insert.
		* WindowWeight(w, x)); //the window value is multiplied in here
}

void SetPseudoCount(Distribution *d)
//gives a reasonable pseudocount
{
	d->pseudocount= MaxWeight(d) / 10000;
}

double GetNormWithWindow(Distribution *dist, int32my gap, Window w)
// get the norm of the prob dist once the window has been applied (and pseudocounts included)
// this is problematic: time overhead here may be high if called for every likelihood estimate
{
	double sum=0;
	int32my i,s;
	s=dist->size;
	// there's a contribution from all values between 0 and len1-len2, but some don;t need doing explicitly
	for (i=0;i < s; ++i) //over all values covered by distribution *leaving out pseudocount*;
		sum += GetWeight(*dist, i-gap, gap,  w, 0); //i here is insert in dist, whereas GetWeight takes the overhang as argument.
	sum += dist->pseudocount * w.len2; //adds the sum of pseudocounts after they have been windowed
	// What is going on here is that the pc can be taken outside the sum, and w.len2 turns out to be the sum of the window values.
	return (sum);
}

double LogLikelihood(Distribution *dist, int32my *samples, int32my numSamples, Window w, int32my gap)
// for an array samples of size numSamples,  find the likelyhood that the values in the array
// came from the distribution dist with window w (which should be set with gap=0).
{
	int32my i,*ptr;
	double ll=0.0;
	double rnorm; // r for reciprical
	rnorm= 1.0 / GetNormWithWindow(dist, gap, w); //unfortunately norm must be calc-ed for every likelihood
	printf("norm= %f \n",rnorm);
	//the samples here represent the length of the end of the scaffolds "between" the reads,
	//so the insert distribution applies to this plus the gap size being tested.
	for (i=0,ptr=samples;i<numSamples;++i,++ptr) {
		ll += log( GetWeight(*dist, *ptr, gap, w, dist->pseudocount) * rnorm );
		//printf("sample %d =%d weight %e \n",i, *ptr, GetWeight(*dist, *ptr, gap, w, dist->pseudocount) * rnorm);
	}
	return(ll);
}

int32my MaxLikelihood(Distribution *dist, int32my *samples, int32my numSamples, Window w, int32my start, int32my end, int32my step)
// for an array samples of size numSamples,  find the max likelihood for the gap size between start and end
{
	int32my i;
	double maxll,ll;
	int32my best;
	maxll= - DBL_MAX;
	SetPseudoCount(dist); //set pseuodcounts to something appropriate
	//printf("find likelihoods \n");
	best=-1;
	for (i=start; i<end; i+=step) {
		ll=LogLikelihood(dist, samples, numSamples, w, i);
		printf("i=%"PRI32" ll %f maxll %f \n",i,ll,maxll);
		if (ll>maxll) {
			best=i;
			maxll=ll;
		}
		//printf("for %d like=%f \n",i,ll);
	}
	return(best);
}

double MultiDistributionLogLikelihood(int32my *samples, int32my numSamples, Distribution **distList, int32my numDists, Distribution **whichDist, Window w, int32my gap)
// for an array samples of size numSamples,  find the likelyhood that the values in the array
// came from the distribution dist with window w (which should be set with gap=0).
{
	int32my i,*ptr;
	double ll=0.0;
	Distribution **dist;

	for (i=0; i<numDists; ++i)
		distList[i] -> tempNorm=1.0 / GetNormWithWindow(distList[i], gap, w); //there should only be 2-5 of these typically
	//printf("norm= %f \n",rnorm[i]);
	//the samples here represent the length of the end of the scaffolds that stick out between the reads,
	//so the gap would be the true insert plus this sample;
	for (i=0,ptr=samples, dist=whichDist;i<numSamples;++i,++ptr,++dist) {
		ll += log( GetWeight( **dist, *ptr, gap, w, (*dist)->pseudocount ) * (*dist)->tempNorm );
		//printf("sample %d =%d weight %e \n",i, *ptr, GetWeight(*dist, *ptr, gap, w, dist->pseudocount) * rnorm);
	}
	return(ll);
}

int32my MultiDistributionMLE(int32my *samples, int32my numSamples, Distribution **distList, int32my numDists, Distribution **whichDist, Window w, int32my start, int32my end, int32my step)
// for an array samples of size numSamples,  find the max likelihood for the gap size between start and end
{
	int32my i;
	double maxll,ll;
	int32my best;
	maxll= - DBL_MAX;
	for (i=0; i<numDists; ++i)
		SetPseudoCount(distList[i]); //set pseuodcounts to something appropriate
	//printf("find likelihoods \n");
	best=-1;
	for (i=start; i<end; i+=step) {
		ll=MultiDistributionLogLikelihood(samples, numSamples, distList, numDists, whichDist, w, i);
		//printf("i=%d ll %f maxll %f \n",i,ll,maxll);
		if (ll>maxll) {
			best=i;
			maxll=ll;
		}
		//printf("for %d like=%f \n",i,ll);
	}
	return(best);
}

int32my LoadInsertSizeSamples(Distribution *dist, char *fname, ContigList *contigList, char cigarTypeChar)
//reads data from a smalt file and finds the distance between pairs located on the same contig.
//smalt output should contain cigar lines with mate pairs sequential.
//Finds the freq distribution of the insert size from this in the Distribution structure dist. dist must be set up with makeEmptyDistribution first.  Samples can be added to existing weight.
// returns number of samples used.
{
        FILE *fp;
	//Contig *contig;
        char line[MaxLineLen],contigName[ReadNameLen],cigar[12],readDir[2]; //halfway houses for file data.
	char cigarType[]="cigar:A";
	char *rval;
        int64my numSamples;
        int32my *weight, insert;
	int32my endPos1,startPos1,endPos2,startPos2;
	int32my pos1,pos2;
	int32my score1,score2;
	int32my index;

	printf("LoadInsertSizeSamples: getting intra-contig pairs to estimate insert distribution\n");
	cigarType[6]=cigarTypeChar; //this can be A B or C.
        //---open file------------------------
        fp = fopen(fname,"r");
        if (fp==NULL) {fprintf(stderr,"Could not open file %s\n",fname); exit(1);}
        //---go through smalt file to find the distance between read pairs
	numSamples=0;
	//first move over preamble 
	//line now contains the first cigar line.
       while (fgets(line,MaxLineLen,fp) != NULL) {
		while(strncmp(line,cigarType,7) && (rval= fgets(line,MaxLineLen,fp))!=NULL );
		if (rval == NULL) break;
                sscanf(line,"%s %*s %*d %*d %s %s %"SCN32" %"SCN32" %*s %*d",
                        cigar,
			readDir,
                        contigName,
			&startPos1,
                        &endPos1);
		//we should find one read going forward, and one backward.  From the forward one we want the end of the match and
		//from the backward one we want the start of the match from which to measure the insert.  Counting from the "head" of the reads
		//will be our standard when estimating the gapsize
		pos1 = *readDir=='+' ? endPos1 : startPos1 ;
                score1=atoi(cigar+8); //first token, "line", should be eg "cigar.D:59" 59 is score
		if (fgets(line,MaxLineLen,fp) == NULL) {fprintf(stderr,"%s file reading %s error: no pair to a cigar:A read.  Is file smalt output? \n",line,fname); exit(1);}
               sscanf(line,"%s %*s %*d %*d %s %s %"SCN32" %"SCN32" %*s %*d",
                        cigar,
			readDir,
                        contigName,
			&startPos2,
                        &endPos2);
		pos2 = *readDir=='+' ? endPos2 : startPos2 ; 
                score2=atoi(cigar+8); //first token, "line", should be eg "cigar.D:59" 59 is score
                if (*(cigar+6)!=cigarTypeChar) {fprintf(stderr,"%s file reading %s error: no pair to a cigar:%c read.  Is file smalt output? \n",line,fname,cigarTypeChar); exit(1);}
		//find out if this read pair should contribute to the insert size weight
                if (score1<MATCH_THRESHOLD || score2<MATCH_THRESHOLD) continue;
		//printf ("contig %s\n",contigName);
		index=GetIndexFromName(contigName, contigList);
		//printf("contig name %s index %" PRI32 "\n", contigName, index);
		if (index==-1) continue;
		//printf ("contig %s len %d \n",contig->name,contig->len);
		// We want the pdf for insert sizes over the genome.  Taking all pairs occuring on short contigs would be biased: more short inserts would counted.
		//contigList->contig[index].len is the length of the contig with index "index".
                if ( (pos1>pos2 ? pos1 : pos2) < dist->size || (pos1>pos2 ? pos2 : pos1) > contigList->contig[index].len - dist->size) continue; 
		++ numSamples; //we now know that we have a valid pair here for our sample.
		if ( (insert= (int32my) abs(pos2 - pos1)) < dist->size )  ++ dist->weight[ insert ];
		//if(abs(pos2 - pos1) < 10000) printf("%d ",pos2 - pos1);
        }
	if (numSamples==0)  printf("------ LoadInsertSizeSamples: no lines in file %s with cigar:%c were used.\n------ Wrong file type, no good pairs or all pairs closer than %"SCN32" bp to ends of contigs.\n",fname,cigarTypeChar,dist->size);
	dist->norm += numSamples;
	printf("# pairs far enough from contig ends to be fair samples: %"SCN64" \n", numSamples );
        //---close file------------------------
        if (fclose(fp)) {fprintf(stderr,"Could not close file %s\n",fname); exit(1);}
	return (numSamples);
}

// ****************************************************
// unimportnant stuff follows which can be removed when not testing.
// some functions for testing and degugging.

double MaxWeightWithWindow(Distribution *d, Window w,int32my gap)
{
	int32my i;
	double weight, *ptr;
	double max = 0.0;
	for (i=0,ptr=d->weight; i<d->size; ++i,++ptr) {
		weight= GetWeight(*d, i - gap, gap, w, 0);
		if (weight > max) max=weight;
	}
	return(max);
}

int32my *GenerateSamples(Distribution *d, int32my numSamples)
//generates samples from the prob dist in the simplest way
//using built-in random numbers.
{
	int32my i, s, *ptr, *samples;
	double max;
	max=MaxWeight(d);
	samples = (int32my *) malloc( sizeof(int32my) * numSamples);
	for (i=0,ptr=samples; i<d->size; ++i,++ptr) {
		do 
			s= rand() % d->size; //a uniformaly random sample
		while ( (rand() / (double) RAND_MAX) * max <= d->weight[s] );
		*ptr=s;
	}
	return(samples);
}

int32my *GenerateSamplesWithWindow(Distribution *d, int32my numSamples, Window w, int32my gap)
//generates samples from the prob dist in the simplest way
//using built-in random numbers.
{
	int32my i, s, *ptr, *samples;
	double max;
	double prob;

	max=MaxWeightWithWindow(d, w, gap);
	d->pseudocount=0.0;

	samples = (int32my *) malloc( sizeof(int32my) * numSamples);
	//printf("generating samples \n");
	for (i=0,ptr=samples; i<numSamples; ++i,++ptr) {
		prob= GetWeight(*d, i - gap, gap, w, 0) ;
		do 
			s= (rand() % d->size) - gap; //a uniformly random sample of *overhangs* (=insert-gap)
		while ( (rand() / (double) RAND_MAX) * max > GetWeight(*d, s, gap, w, 0) );
		*ptr=s;
	}
	return(samples);
}

void PrintWeightsWithWindow(Distribution *d, Window w, int32my gap)
{
	int32my i, *ptr;
	double rnorm; // r for reciprical
	rnorm= 1 / GetNormWithWindow(d, gap, w);
	printf("Distribution for gap size %"PRI32" \n",gap);
	for (i=0;i < d->size; ++i)
		printf(" overhang %"PRI32" : weight %e win %e\n",i,((i+gap)>=0 && (i+gap)<d->size) ? d->weight[i + gap] + 0 : 0,GetWeight(*d, i, gap,  w, d->pseudocount) * rnorm);
}

/*************************************************************************************************************************/
////////////////////////////
