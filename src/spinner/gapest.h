// gapest.c
// routines for finding the distribution of insert sizes
// and using these to estimate gap sizes

struct _ContigList;

typedef struct
{
	int32my len1, len2;
} Window;

typedef struct _Distribution//to store a frequency histogram or pdf
{
        double *weight;    // weight of this value (e.g. prob or frequency)
        double norm; // normalisation (e,g, 1.0 or the number of samples)
	int32my size;       // the size the the samples array
	double pseudocount; //pseudocount can be added to weights.
	double tempNorm; // ignorethis
} Distribution;

typedef struct //to store a frequency histogram or pdf
{
        Distribution *dist; //array of distributions
	int32my numDist;
} DistSet;

typedef struct //to store a frequency histogram or pdf
{
        int32my **samples; //array of pointers to samples for each dist
	int32my *numSamples; //size of the arrays
	int32my numDist;
} SampleSet;

Window GetWindow(int32my len1, int32my len2);

void PrintWeightsWithWindow(Distribution *d, Window w, int32my gap);

double MeanValue(Distribution *d);
//gets the mean from the distribution of samples

double StdDev(Distribution *d);

double MaxWeight(Distribution *d);

int32my *GenerateSamples(Distribution *d, int32my numSamples);

int32my *GenerateSamplesWithWindow(Distribution *d, int32my numSamples, Window w, int32my gap);

Distribution GetEmptyDistribution(int32my size);
//set a distribution of size size with no samples recorded

void DeleteDistribution(Distribution dist);
//free memory for dist

void Normalise(Distribution *dist);
//Makes a probability distribution by normalising dist to 1

void ResetNorm(Distribution *dist);
//The norm should the the sum of the weights
//This forces that to be true (corrects floating point errors etc.)

void Smooth(Distribution dist, double minBlockWeight);
// intended to get rid of zeros of distributions.
// blocks of values in the distribution with aggregate weight less than minBlockWeight are given the mean weight for the block

int32my LoadInsertSizeSamples(Distribution *dist, char *fname, struct _ContigList *contigList, char cigarTypeChar);
//reads data from a smalt file and fins the distance between pairs located on the same contig.  Finds the distribution  of the insert size from this in the Distribution structure dist. dist must be set up with makeEmptyDistribution first.  Samples can be added to existing samples.

int32my MaxLikelihood(Distribution *dist, int32my *samples, int32my numSamples, Window w, int32my start, int32my end, int32my step);
// for an array samples of size numSamples,  find the max likelihood for the gap size between start and end

void SetNormalDistribution(Distribution *dist, double mean, double sigma);

int32my MultiDistributionMLE(int32my *samples, int32my numSamples, Distribution **distList, int32my numDists, Distribution **whichDist, Window w, int32my start, int32my end, int32my step);