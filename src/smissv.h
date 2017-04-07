
#ifndef SMISSVHEADER

#define SMISSVHEADER

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <cstring>
#include <string>
#include <ctype.h>
#include "fasta.h"
#include <vector>
#include <iostream>
#include <algorithm> 
#include <sstream>



#define Max_N_NameBase 400 


static int yes=1;
static int no=0;
static int pri=1;


int sign(int num)
{
  if(num>=0) return 1;
  else return -1;
}


template<class Iterator>
Iterator Unique(Iterator first, Iterator last) 
{ // eliminate repetitive elements from sorted vector
  while (first != last)
    {
      Iterator next(first);
      last = std::remove(++next, last, *first);
      first = next;
    }
  return last;
}


void Print (const std::vector<std::string> &v){
  
  for (int i=0; i<v.size();i++)
    std::cout << i << " " << v[i]<< std::endl;
  std::cout << " Vector size " << v.size() << std::endl;
}


void Printi (const std::vector<int> &v){
  
  for (int i=0; i<v.size();i++)
    std::cout << i << " " << v[i]<< std::endl;
  std::cout << " Vector size " << v.size() << std::endl;
}

char * strrev(char *str)
{
  int i = strlen(str)-1,j=0;
  char ch;
  while(i>j)
    {
      ch = str[i];
      str[i]= str[j];
      str[j] = ch;
      i--;
      j++;
    }
  return str;
}


// *********************************//
int pyplot(char *nameout, char *pdfname)
// *********************************//
{ 
  char syscmd[2000];
  memset(syscmd,'\0',2000);
  sprintf(syscmd,"python ~/bin/software/smissv/src/plot.py %s %s",
	  nameout,pdfname);
  if(system(syscmd) == -1)
    printf("System command error:\n");

  return 0;
}


/* creat char matrix with subscript range cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
	  {
	    printf("error cmatrix: calloc error No. 1 \n");
	    return(NULL);
	  }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
	  {
	    printf("error cmatrix: calloc error No. 2 \n");
	    return(NULL);
	  }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
	  cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}


// numeric to string
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}
// string of numbers to int 
template <class T1>
inline int to_int (const T1& t)
{
  int ii;
  std::stringstream ss(t);
  ss >> ii;
  return ii;
}

struct IdxCompare
{
    const std::vector<int>& target;

    IdxCompare(const std::vector<int>& target): target(target) {}

    bool operator()(int a, int b) const { return target[a] < target[b]; }
};





#endif
