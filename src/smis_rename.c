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
/****************************************************************************/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "smissv.h"



static B64_long sBase;

/* SSAS default parameters   */

static int num_reads=0;
static int n_contigs = 0;

fasta *seq;


int main(int argc, char **argv)
{
    FILE *namef;
    int i,nSeq,args;
    char line[2000]={0};
    fasta *seq;  // double declaration???
  
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Rename_Process(char **argv,int args,int nRead);

    seq=NULL;

    if(argc < 2)
    {
      printf(" %s <input fasta/q> <output fasta >\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
 

    num_reads = 0;
    Rename_Process(argv,args,num_reads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   Subroutine to process alignment information */
/* ====================================================  */
void Rename_Process(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,rc,nSeq = nRead;
     char **DBname,*ptr,RC,nametag1[100],nametag2[100];
     char *line,*st,*ed;
     FILE *fp,*namef,*namef2,*fpOutfast,*fpOutfast2;
     B64_long *read_offsets,big_num,sum_bases;
     int n_patch,idd,stopflag,num_gcs,*contig_index;
     float rate;
     fasta *segg,*seqp;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((seq=decodeFastq(argv[args],&num_seqque,&sBase,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     n_contigs = num_seqque;
     fastaUC(seq,n_contigs);

     nSeq = n_contigs;
     nRead = n_contigs; 

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
         printf("ERROR main:: reads group file \n");
         exit(1);
     }
/*   output all the fastq files    */
     for(i=0;i<nSeq;i++)
     {
        char *dpp;
        int gen_len,nline;

        seqp=seq+i;
        gen_len = seqp->length;
        fprintf(namef,">contig%dsize%d\n",i+1,gen_len);
        //fprintf(namef,">contig/%d/0_%d RQ=0.8\n",i+1,gen_len);

        //fprintf(namef,">%d\n",i+1);

        nline = gen_len/60;
        for(k=0;k<nline;k++)
        {
           for(j=0;j<60;j++)
              fprintf(namef,"%c",seqp->data[k*60+j]);
           fprintf(namef,"\n");
        }
        for(j=0;j<(gen_len-(nline*60));j++)
           fprintf(namef,"%c",seqp->data[nline*60+j]);
        if((seqp->length%60)!=0)
          fprintf(namef,"\n");
     }
     fclose(namef);
}



