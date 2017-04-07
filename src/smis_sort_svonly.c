/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2013, 2014  Genome Research Ltd.                          *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of smis_pipeline.                                     *
 *                                                                          *
 *  ssaha_pileup is free software: you can redistribute it and/or modify it *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/

#include "smissv.h"
#include <api/BamReader.h>
#include <api/SamHeader.h>
#include <api/SamSequenceDictionary.h>
#include <api/BamAux.h>
#include <ctime>

using namespace BamTools;

struct ALIGNMENT {
  std::string refname,readname,nmscore;
  int32_t ontmask,hitrcdex,hitlocus1,hitlocus2,readl,readpos;
  long int ontrdlen,ontrdidt,ontrdtag;
  int superl;
  uint32_t hitread1,hitread2;
  uint16_t mapscore;
  float identy;

  bool operator() (ALIGNMENT i, ALIGNMENT j) { 
    return( ((i.ontrdlen<<46) + (i.ontrdtag<<32) + (i.ontrdidt<<8) +  i.readpos)  
    	    <  ((j.ontrdlen<<46) + (j.ontrdtag<<32) + (j.ontrdidt<<8) + j.readpos)); 
  }
} alignment;

static  std::vector<ALIGNMENT> myaligns;


int GetIntTag(const BamAlignment &aln, const std::string &tag, unsigned long int & mytag);
int ReadBwa(std::string infile);
int ReadSmalt(FILE *namef);

int Sort_and_Select();


static char nameout[100];
static  FILE *outfile;
//static  FILE *namef;

static int Ncontigs;
static int Fcontigs;
static  std::vector<int> finalctg;

static int set_score = 10; // mapping score
static float set_identy=60.0;  // id 
static  std::vector<std::string> listctg;
static  int nseq=0;
static int printsome=5;



int main(int argc, char **argv)
{
  int args;
  int readok=-1;
  FILE *namef;
 
  time_t starttime = time(0);
  std::cout << "Starting at " << starttime <<  std::endl;

  if(argc < 2){
    printf("Usage: %s [-identy 60.0] [-score 10] <alignment_file>\n",argv[0]);
    exit(1);
  }

 args=1;
 for(int i=1;i<argc;i++){
   if(!strcmp(argv[i],"-score")){
     sscanf(argv[++i],"%d",&set_score);
     args=args+2;
   }else if(!strcmp(argv[i],"-identy")){
     sscanf(argv[++i],"%f",&set_identy);
     args=args+2;
   }
 }
  
 myaligns.reserve(100000);


 std::string infile=argv[args];
 bool smalt=1;
 std::string ftype=infile.substr(infile.size()-3);
 if(ftype=="bam")
   smalt=0;
 
 
 if(smalt){
    namef = fopen(argv[args],"r");
    readok=ReadSmalt(namef);
 }else{
   readok=ReadBwa(infile);
 }
 
 time_t readdone = time(0);
 std::cout << "Done reading at " << readdone
	   << " Total time so far " << readdone-starttime<<  std::endl;
 
 if(readok == 0){
   memset(nameout,'\0',100);          
   sprintf(nameout,"mates_id%d_mscore%d.out",
	   set_identy,set_score);
   Sort_and_Select();
 }else{
   printf("\nERROR:: Alignment file not properly read!! \n");
   exit(1);
 }

 printf("Job finished for %d reads from file %s\n",nseq,argv[args]);
 
 time_t finish = time(0);
 std::cout << "All done at " << finish
	    << " Second part took " << finish-readdone 
	    << " sec \nTotal time " << finish-starttime << " sec"<<  std::endl;


 return(0);
  
} /* end of the main */

  
// ******************* //
int Sort_and_Select()
// ******************* //
{
  int i1,i2;

  if((outfile = fopen(nameout,"w")) == NULL){
    printf("\nERROR:: cannot write on %s \n",nameout);
    exit(1);
  }

  // sort
  std::sort(myaligns.begin(), myaligns.end(), alignment);

  listctg.erase( Unique( listctg.begin(), listctg.end() ), listctg.end() );
  Ncontigs= listctg.size();
  printf("Total number of contigs: %d \n", Ncontigs);


  for (int j=0; j<myaligns.size();j++){

   if(j%2 == 0){
     i1 = j;
   }else{
     i2 =j;
     
     
     // to avoid cases where not reading the whole alignemnt file
     if( myaligns[i1].readname.substr(0,myaligns[i1].readname.length()-2) 
	 !=myaligns[i2].readname.substr(0,myaligns[i2].readname.length()-2)){
       printf("ERROR:: read mate not paired? %d %d %s %s\n",i1,i2,
	      myaligns[i1].readname.c_str(),
	      myaligns[i2].readname.c_str() );
       exit(1);
      }
     
      // difference from smis_pipeline: was map_score[i1] > set_score instead of >=
     if( (myaligns[i1].mapscore >= set_score) && (myaligns[i2].mapscore >= set_score)
	 && (myaligns[i1].identy >=set_identy) && ( myaligns[i2].identy >= set_identy)
	 && (myaligns[i2].ontmask > 0)) {
	
	
 	std::string sc = myaligns[i1].nmscore;
	if(myaligns[i2].mapscore < myaligns[i1].mapscore)
	  sc=myaligns[i2].nmscore;

       
       int distance=0;
       if((myaligns[i1].hitrcdex == 0)&&(myaligns[i2].hitrcdex == 1))
	 distance=myaligns[i2].hitlocus2-myaligns[i1].hitlocus1-myaligns[i2].ontrdlen;
       else if((myaligns[i1].hitrcdex == 1)&&(myaligns[i2].hitrcdex == 0))
	 distance=myaligns[i1].hitlocus2-myaligns[i2].hitlocus1-myaligns[i2].ontrdlen;
      
       
 	if(yes)fprintf(outfile,"cigar:%s %s %s %d %d %d %d %d %5.2f %d %d || %s %5.2f %d %d %d\n",
		       sc.c_str(),myaligns[i2].refname.c_str(), myaligns[i1].refname.c_str(),
		       myaligns[i1].hitread1,  myaligns[i1].hitread2,
		       myaligns[i1].hitlocus1, myaligns[i1].hitlocus2,
		       myaligns[i1].hitrcdex,myaligns[i1].identy,
		       myaligns[i1].superl, myaligns[i1].ontrdlen,
		       myaligns[i2].readname.c_str(),
		       myaligns[i2].identy,
		       myaligns[i2].hitlocus1,myaligns[i2].hitlocus2,
		       distance);
	
	
       // eliminate/reduce this next part??
	  
       std::string chr1=myaligns[i2].refname;
       chr1.erase(chr1.find("size"));
       
       int ichr1;
       std::string sub;
       if(chr1.size()==7)
	 sub = chr1.substr(chr1.size()-1, 50);
       else if (chr1.size()==8)
	 sub = chr1.substr(chr1.size()-2, 50);
       else if (chr1.size()==9)
	 sub = chr1.substr(chr1.size()-3, 50);       
       ichr1=std::atoi(sub.c_str());
       

       finalctg.push_back(ichr1);
       
       
     }//selection
   }    
 }
 
 std::sort (finalctg.begin(),  finalctg.end());
 finalctg.erase(Unique(finalctg.begin(), finalctg.end() ),finalctg.end() );
 
 Fcontigs= finalctg.size();
 
 
 fprintf(outfile,"%d  %d ",Ncontigs,Fcontigs);
 for(int k=0;k<Fcontigs;k++){    
   fprintf(outfile," %d ",finalctg[k]);
 }
 fprintf(outfile,"\n");
 fclose(outfile);


 return(0);
}



// ******************* //
int ReadSmalt(FILE *namef)
// ******************* //
{
  int read=1;
  std::string prevchr;
  nseq=-1;

  while(read){

    if(nseq<printsome)pri=1;
    else pri=0;

    char nmscore[250],rname[250];
    char tempc1[60],ctgname[60],nont_RC[60];
    int  nhit_read1,nhit_read2,nhit_locus1,nhit_locus2;
    int allength,superl;
    float nidenty;
    
    ALIGNMENT thisal;
    nseq++;
   if(fscanf(namef,"%s %s %s %s %d %d %d %d %s %d %f %s %d %d",
	      nmscore,tempc1,rname,ctgname,  
	      &nhit_read1,&nhit_read2,
	      &nhit_locus1,&nhit_locus2,
	      &nont_RC,&allength,&nidenty,
	      tempc1,&superl)!=EOF){
      
 
      std::string thischr=ctgname;
      if(prevchr!=thischr){
	listctg.push_back(thischr);
	prevchr=thischr;
      }

      
      // read 1 or 2
      char *st;
      thisal.readname = rname;
      thisal.readpos = to_int(thisal.readname.substr (thisal.readname.find_last_of("/")+1));
      //thisal.readname = thisal.readname.substr(0,thisal.readname.length()-2);


      // Alignment Mapping score
      thisal.nmscore=nmscore; 
      thisal.nmscore=thisal.nmscore.erase(0,10);
      std::string tempstr=thisal.nmscore.substr(0,thisal.nmscore.find_last_of(":"));
      thisal.mapscore = to_int(thisal.nmscore.substr (thisal.nmscore.find_last_of(":")+1));
     
      thisal.ontmask = 0;
      if(tempstr == "A")
	thisal.ontmask = 1;
      else if(tempstr == "D")  // reads mapped in different chr or contigs, useful for scaffolding
	thisal.ontmask = 2;
     
      thisal.hitlocus1 = nhit_locus1; // leftmost pos in ref
      thisal.hitlocus2 = nhit_locus2;
      thisal.hitread1 = nhit_read1; // leftmost pos in ref
      thisal.hitread2 = nhit_read2;
      thisal.ontrdtag = to_int(thisal.readname.substr(10,6));	//insert length	  
      thisal.ontrdlen = to_int(thisal.readname.substr(26,7));
      thisal.ontrdidt = to_int(thisal.readname.substr(17,8));

      thisal.superl = superl;
      thisal.identy = nidenty;
      thisal.refname = ctgname;
      if(nont_RC[0] == 'F'){
	thisal.hitrcdex = 0;
      } else {
	thisal.hitrcdex = 1;
      }

      if(pri)std::cout <<  nmscore << " " << thisal.readname << " " << thisal.readname.substr(10,6) << " " 
		<< thisal.readname.substr(17,8) << std::endl;

      
      //if(nseq>10)read=0;
      myaligns.push_back(thisal);
    }
    thisal = ALIGNMENT(); // re-initialize thisal
    if(feof(namef))read=0;  
  }
  fclose(namef);


  return(0);
}




// *************************************************************************************** //
int ReadBwa(std::string infile)
// *************************************************************************************** //
 {
  BamReader reader;
   if ( !reader.Open(infile) ) {
     std::cerr << "Could not open input BAM files." << std::endl;
     exit(1);
   }
 
   const SamHeader header = reader.GetHeader();
   const RefVector references = reader.GetReferenceData();

  BamAlignment al;
  int read=1;
  std::string prevchr;
  nseq=-1;
 
  while (read &&  reader.GetNextAlignment(al)) {
    
    if(nseq<printsome)pri=1;
    else pri=0;


    if(al.IsPrimaryAlignment()  && al.IsMapped() 
       && al.IsMateMapped() && al.IsProperPair()) {
      //&& al.RefID == al.MateRefID){

      nseq++;

      std::string thischr=to_string(al.RefID);
      if(prevchr!=thischr){
	listctg.push_back(thischr);
	prevchr=thischr;
      }
          
      ALIGNMENT thisal;

      if(al.IsReverseStrand()){
	thisal.hitrcdex = 1;
      } else {
	thisal.hitrcdex = 0;
      }
      thisal.mapscore = al.MapQuality;
      thisal.readname = al.Name + "/?";
      thisal.hitlocus1 = al.Position; // leftmost pos in ref
      thisal.hitlocus2 = al.GetEndPosition();  // rightmost pos in ref ? or in read??? check
      thisal.refname = references[al.RefID].RefName;
      thisal.ontrdtag = to_int(al.Name.substr(10,6));		  
      thisal.ontrdlen = to_int(al.Name.substr(26,7));
      thisal.ontrdidt = to_int(al.Name.substr(17,8));
      thisal.readl = al.Length;
      
    
      std::string strtag = "NM";
      unsigned long int value;
      if(GetIntTag(al,strtag,value)){
    	printf("WARNING:: cannot find \"NM\" field in your bam file, each read identity will be set to 99\%  \n");
	value=al.AlignedBases.size()-1;
      }
      thisal.identy = (value*100.) /(al.AlignedBases.size());
      
      //ctgname no need
      //ont_RC[i] no need
      //superlength[i])  not needed?
       
      // not defined in bwa
      thisal.hitread1 = 0; // leftmost pos in read
      thisal.hitread2 = 0; // rigthmost pos in read
      thisal.ontmask = 1;  // A cut: only properly paired reads
      thisal.readpos = 0;  // no way to distinguish read 1 and 2 from bwa output
      thisal.nmscore= "A:"+ to_string(thisal.mapscore);   
      thisal.superl=0;  //superlength[i])  not needed?
   

      myaligns.push_back(thisal);
      thisal = ALIGNMENT(); // re-initialize thisal
     

    }//selection
    //if(nseq>1000)read=0;
  }//while


  return(0);
}
// *************************************************************************************** //
int GetIntTag(const BamAlignment &aln, const std::string &tag, unsigned long int & mytag){
// *************************************************************************************** //
  char type;
  uint32_t val;

  if (!aln.HasTag(tag)){
    std::cout << "ERROR:: no tag found " << std::endl;
    return(1);
  }
  if (!aln.GetTagType(tag,type)){
    std::cout << "ERROR:: tag " << tag << " not found " << std::endl;
    return(2);
  }
  if(aln.GetTag(tag,val)){
    mytag = val;
    return(0);
  }else{
   std::cout << "ERROR:: tag " << tag << " has wrong format " << std::endl;
   return(3);
  }
}
