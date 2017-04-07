//constants and such like needed in many places.

#include <inttypes.h>

//portable data types
//NB: I assume in most places that no scaffold will
//contain more than 2147483648 bases, so things to do with placement are stored in 32 bit ints.

#define SCN64 SCNdLEAST64
#define SCN32 SCNdLEAST32
#define PRI64 PRIdLEAST64
#define PRI32 PRIdLEAST32
#define PRI16F PRIdFAST16
#define PRIBOOL PRIdLEAST8
#define PRIBOOLF PRIdFAST8


#define uint16my uint_least16_t
#define uint32my uint_least32_t
#define uint64my uint_least64_t
#define int16my int_least16_t
#define int32my int_least32_t
#define int64my int_least64_t
#define boolmy int_least8_t
#define boolmy_fast int_fast8_t

#define uint16f uint_fast16_t
#define uint32f uint_fast32_t
#define uint64f uint_fast64_t
#define int16f int_fast16_t
#define int32f int_fast32_t
#define int64f int_fast64_t


 #define STORE_READ_NAMES
// this includes the readname data in the data loaded into memory.

#define DEFAULT_FILES_NAME "files.txt"
#define DEFAULT_SETTINGS_NAME "settings.txt"
#define DEFAULT_FASTQ_OUT_NAME "(NOT SAVING IN THIS FORMAT)"
#define NO_NAME "(NO SAVE)"


#define INSERT_SIZE_TO_STD_DEV_RATIO 5
//the std dev assigned to reads insert size is the given insert size divided by this at the moment.
#define DEFAULT_REPLACEMENT_SIZE_FOR_NEGATIVE_GAP 10
//used when writing out the contigs in fastq / fasta format at the moment.
#define FRAYED_END_SIZE 1500
//should get rid of this rubbish really
#define MIN_INCERTAINTY_REC 1.0
//the standard deviation of an edge is not allowed to be any set smaller than 1/this.
#define ReadNameLen 60
//ReadNames can't be longer than this.
#define MaxLineLen 1024
//when reading files


// one day the following will be needed to align ends.
#define MAX_SEQ_LENGTH 6000
//max size for the part of contigs to be aligned
#define NO_SKIP -1000000
// special value for skipLen find closest parameter indictating no skip test.
