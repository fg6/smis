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

/* SSAS default parameters   */
static int   set_score = 10; // mapping score
static float set_identy=60.0;  // id 
static int   set_inslen = 2500; // mapping score
static int selcontig = 1; 
static float ratio = 0.5;
   
static std::vector<int> sorted_index;
static  std::vector<float> vnormins;


static FILE *namef; 
static FILE *outfile;
static char nameout[100];

static char fortest[100];
static FILE *outtest1;
static FILE *outtest2;

const int debug=1;

int read_and_print(void);
int norm_ont(void);
int printtest(void);

static  std::string chromosome="contig";
//static  std::vector<int> vmapsc,vid1,vid2;
static  std::vector<int> vinslen,vposition;
static  std::vector<std::string> vsmreads;




// *********************************//
int main(int argc, char **argv)
// *********************************//
{
  int args;
  if(argc < 2){
    printf("Usage: %s [-identy 60.0] [-score 10] [-chr 1] [-ilen 2500] <alignment_file>\n",
	   argv[0]);
    exit(1);
  }
  args=1;
  for(int i=1;i<argc;i++)
    {
      if(!strcmp(argv[i],"-ilen")){
	sscanf(argv[++i],"%d",&set_inslen); 
	args=args+2;
      }else if(!strcmp(argv[i],"-score")){
	sscanf(argv[++i],"%d",&set_score);
	  args=args+2;
      }else if(!strcmp(argv[i],"-identy")){
	sscanf(argv[++i],"%f",&set_identy);
	args=args+2;
      }else if(!strcmp(argv[i],"-chr")){
	sscanf(argv[++i],"%d",&selcontig);
	  args=args+2;
      }else if(!strcmp(argv[i],"-ratio")){
	sscanf(argv[++i],"%f",&ratio);
	  args=args+2;
      }
    }
  


  if((namef = fopen(argv[args],"r")) == NULL){
    printf("ERROR main:: missing alignment file !! \n");
      exit(1);
  }

  std::stringstream mystring;
  mystring << selcontig;
  chromosome += mystring.str();
  if(yes) std::cout << "\n Analyzing Chromosome " << chromosome 
		    << " with min identity score " << set_identy
		    << ", min mapping score " << set_score
		    << ", cut on difference of insert length " << set_inslen
	    	    << std::endl;
  

  if(debug) printtest();
  read_and_print(); 
  
 
  norm_ont();
  if(debug) fclose(outtest1);
  if(debug) fclose(outtest2);


  //Write output
  memset(nameout,'\0',100);          
  sprintf(nameout,"chr%d.dat",selcontig);
  if((outfile = fopen(nameout,"w")) == NULL){
    printf("ERROR:: cannot write on %s \n",nameout);
    exit(1);
    }
  for(int v=0; v<vposition.size(); v++){
    int ix=sorted_index[v];
    fprintf(outfile,"%s %d %d %f\n",vsmreads[ix].c_str(),
	    vposition[ix],vinslen[ix],
	    sign(vinslen[ix])*vnormins[v]);
  }
  fclose(outfile);
  

  char pdfname[100];
  sprintf(pdfname,"chr%d",selcontig);
  pyplot(nameout,pdfname);


  if(yes)printf(" Total number of SV found: %d, listed in file %s plotted in pdf %s.pdf \n", 
		vposition.size(),nameout,pdfname); 

  
  return EXIT_SUCCESS;

} // end of main


// *********************************//
int norm_ont()
// *********************************//
{ 

  vnormins.push_back(-10);
  int ix=sorted_index[0];
  if(no)fprintf(outtest2,"%d %d %d %f\n",0,vposition[ix],vinslen[ix],
		   vnormins[0]);

  for(int i=1; i<vinslen.size();i++){
    int ix=sorted_index[i];
    int ixm1=sorted_index[i-1];

    float idd=1;
    float rate,rate2;
      
    if((vposition[ix] - vposition[ixm1]) != 0)
      idd = vposition[ix] - vposition[ixm1];
   
    rate = fabs(idd/vinslen[ix]);
   
    if((1./rate) <= 1.0)   // eliminate positions for which distance in genome > insert size
      rate = 1.0;
    rate2 = 1.0-pow(rate,0.3);

    if(rate2 >= ratio){
      vnormins.push_back(rate2);
      if(debug)fprintf(outtest2,"%d_%d_%1.0f_%f %d %f\n",vinslen[ix],vposition[ix],
		     idd,1/rate,vposition[ix],sign(vinslen[ix])*vnormins[i]);
      if(no)fprintf(outtest2,"%d %d %d %f\n",i,vposition[ix],vinslen[ix],
		       vnormins[i]);
    }else{
      vnormins.push_back(-10);
      if(no)fprintf(outtest2,"%d %d %d %f\n",i,vposition[ix],vinslen[ix],
		       vnormins[i]);
    }
    
  }

}

// *********************************//
int read_and_print()
// *********************************//
{ 
  int nseq=-1;
  int nchr=-1;
  int npassed=-1;
  int read=1;
  char split_char= ' ' ;
  char line[2000]={0};

  pri=0;
  int is=0;
  int fs=0;

  fgets(line,2000,namef);
 


  while(read) 
    {
      nseq++;
      if(nseq>is && nseq<fs)pri=1;
      else pri=0;      

      

      fgets(line,2000,namef);
           
      if(no) 
	printf("line \n %d %s\n",nseq,line);
      
      std::vector<std::string> words;
      std::istringstream split(line);
      
      for(std::string each; getline(split, each, split_char); 
	  words.push_back(each));
           
      if(pri) {
	for(int ii=0; ii<words.size();ii++)
	  std::cout << " splittato " << ii << " " <<  words[ii] << " " << std::endl;
	std::cout << std::endl;
      }

      std::string thischr1=words[1] ;
      std::string thischr2=words[2] ;
      thischr1.erase(thischr1.find("size"));
      thischr2.erase(thischr2.find("size"));



      std::string ms = words[0].substr (words[0].find_last_of(":")+1);
      int mapsc = std::atoi(ms.c_str());
     
      std::string mapq = words[0].substr (words[0].find(":")+1);
      mapq.erase(mapq.find(ms));
      mapq.erase(1,1);


      int id1=std::atoi(words[8].c_str());
      int id2=std::atoi(words[13].c_str());
      int inslen=std::atoi(words[16].c_str());
      int position=std::atoi(words[5].c_str());


      const char *smreads= words[12].c_str();
      const char *mmq=mapq.c_str();
	


      if( (chromosome==thischr1  || chromosome==thischr2))
	nchr++;

      if( (chromosome==thischr1  || chromosome==thischr2) 
	  && fabs(inslen) >=  set_inslen
	  && id1 >= set_identy &&  id2 >= set_identy
	  && mapsc >= set_score
	  && mapq == "A"
	  ){
	
	
	if(npassed<10 && no) 
	  std::cout << " found chr! " << chromosome
		    << thischr1 << " " << thischr2 
		    << " " << id1<< " " << id2 << std::endl;

	//vmapsc.push_back(mapsc);
	//vid1.push_back(id1);
	//vid2.push_back(id2);
	vinslen.push_back(inslen);
	vposition.push_back(position);
	vsmreads.push_back(words[12]);
	

	
	npassed++;
      }
    
      if(feof(namef))read=0;  
    }
  fclose(namef); 
  pri=1;

  if(nchr==-1){
    std::cout << "\n  Error !! Chromosome/contig " 
	      << chromosome << " not found!" << std::endl;
    return 1;
  }else if(npassed==-1){
    std::cout << "\n  !!! Error !!  Chromosome/contig " 
	      << chromosome << " found, but no reads passed the criteria:\n   "
	      << " min identity score " << set_identy
	      << ", min mapping score " << set_score
	      << ", cut on difference of insert length " << set_inslen
	      <<  std::endl;
    return 1;
  }

  
  


  // Sort index
  for(int jj=0; jj<vposition.size();jj++)
    sorted_index.push_back(jj);
  std::sort(sorted_index.begin(),sorted_index.end(),
	    IdxCompare(vposition));



  if(debug){
    for(int v=0; v<vposition.size(); v++){
      int ix=sorted_index[v];
      fprintf(outtest1,"%s %d %d\n", vsmreads[ix].c_str(),
	      vposition[ix],vinslen[ix]);
    }
  }//debug
  


  return 0;
  
} /* end of the r_and_p */



// ********************************* //
int printtest()
// ********************************* //
{ 
  memset(fortest,'\0',100);   
  sprintf(fortest,"chr%d-indels.dat",selcontig);
  if((outtest1 = fopen(fortest,"w")) == NULL){
    printf("ERROR:: cannot write on %s \n",fortest);
    exit(1);
  }
  
  memset(fortest,'\0',100);   
  sprintf(fortest,"chr%d-indels-%2.1f.dat",selcontig,ratio);
  if((outtest2 = fopen(fortest,"w")) == NULL){
    printf("ERROR:: cannot write on %s \n",fortest);
    exit(1);
  }
 

}

