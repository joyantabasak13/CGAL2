#include <iostream>
#include <string>
#include <limits>
#include <vector>
using namespace std;

#include "cgal.h"
#include <math.h>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 500
#define MAX_NAMELENGTH 500

struct SAM
{
	char qname[MAX_NAMELENGTH];
	int flag;
	char rname[MAX_NAMELENGTH];
	int pos;
	int mapq;
	char cigar[MAX_READLENGTH];
	char rnext[MAX_NAMELENGTH];
	int pnext;
	int tlen;
	char seq[MAX_NAMELENGTH];
	char qual[MAX_NAMELENGTH];
	char md[MAX_NAMELENGTH];
	int ih;
	int readLength;
	int ncount;
};


FILE *outFile;
FILE *mapFile;
FILE *unFile;
FILE *statFile;

vector <SAM *> reads1;
vector <SAM *> reads2;

char line1[MAX_REC_LEN];
char line2[MAX_REC_LEN];

long int unCount=0;
long int totalCount=0;
int maxReadLength=0;
double insertSizeMean=0;
double insertSizeVar=0;
double squaredError=0;

long int MAX_FRAGMENT_SIZE=5000;

void reverse(char *reverse, char *read)
{

	char ch='A';	
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else		
			reverse[readLength-i]='N';						
	}
	reverse[readLength]='\0';

}




void writeSam(SAM* read, FILE *out)
{
/*	fprintf(out,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->rname,read->pos,
		read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
*/

	fprintf(out,"%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
		read->cigar,read->tlen,read->seq,read->md,read->ih);

}

void printSam(SAM* read)
{
	printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->rname,read->pos,
		read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);


/*	printf("%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
		read->cigar,read->tlen,read->seq,read->md,read->ih);
*/
}

char temp[MAX_READLENGTH];

SAM *read1, *read2;
int readLength1,readLength2;
int pos1, pos2;
int insertSize;

void printVectors(FILE * out, FILE * u)
{
		
	
	int ih=reads1.size();

	for(int i=0;i<ih;i++)
	{
		read1=reads1[i];
		read2=reads2[i];

		
		if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0 || strcmp(read1->rname,read2->rname)!=0)
		{
		
			
			if(read1->readLength>maxReadLength)
				maxReadLength=read1->readLength;
	
			if(read2->readLength>maxReadLength)
				maxReadLength=read2->readLength;
	

			if(read1->ncount/(double)read1->readLength < 0.8 && read2->ncount/(double)read2->readLength < 0.8)
			{
				unCount++;
				totalCount++;

			
				fputs("@",u);
				fputs(read1->qname,u);
				fputs("\n",u);
			

				int strandNo=(read1->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read1->seq);
					fputs(temp,u);
				}
				else	
				{				
					fputs(read1->seq,u);
				}
				fputs("\n",u);

				fputs("+",u);
				fputs(read1->qname,u);
				fputs("\n",u);

				fputs(read1->qual,u);
				fputs("\n",u);

				fputs("@",u);
				fputs(read2->qname,u);
				fputs("\n",u);
			
				strandNo=(read2->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read2->seq);
					fputs(temp,u);
				}
				else
				{				
					fputs(read2->seq,u);
				
				}
				fputs("\n",u);

				fputs("+",u);
				fputs(read2->qname,u);
				fputs("\n",u);

				fputs(read2->qual,u);
				fputs("\n",u);
				
				for(int i=0;i<reads1.size();i++)
					delete reads1[i];
			
				for(int i=0;i<reads2.size();i++)
					delete reads2[i];
	

				reads1.clear();
				reads2.clear();
				
				return;
			}

		}
		else
		{

			pos1=read1->pos;
			if(read1->readLength>maxReadLength)
				maxReadLength=read1->readLength;
			read1->ih=ih;
	
	
			pos2=read2->pos;
			if(read2->readLength>maxReadLength)
				maxReadLength=read2->readLength;
	
			read2->ih=ih;
	
		//	writeSam(read1,out);
		//	writeSam(read2,out);
						
			SAM *read=read1;
			fprintf(out,"%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
		read->cigar,read->tlen,read->seq,read->md,read->ih);

			read=read2;
			fprintf(out,"%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
		read->cigar,read->tlen,read->seq,read->md,read->ih);


		}
	}
	totalCount++;

	
			
	for(int i=0;i<reads1.size();i++)
		delete reads1[i];

	for(int i=0;i<reads2.size();i++)
		delete reads2[i];
	

	reads1.clear();
	reads2.clear();

}



SAM *getSAM(char *line)
{
	SAM *sam=new SAM;

	
	int i,j;

	i=0;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->qname[j]='\0';
			break;			
		}
		sam->qname[j++]=line[i++];
	}

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			temp[j]='\0';
			break;			
		}
		temp[j++]=line[i++];
	}
	sam->flag=atoi(temp);

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->rname[j]='\0';
			break;			
		}
		sam->rname[j++]=line[i++];
	}

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			temp[j]='\0';
			break;			
		}
		temp[j++]=line[i++];
	}
	sam->pos=atoi(temp);

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			temp[j]='\0';
			break;			
		}
		temp[j++]=line[i++];
	}
	sam->mapq=atoi(temp);

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->cigar[j]='\0';
			break;			
		}
		sam->cigar[j++]=line[i++];
	}

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->rnext[j]='\0';
			break;			
		}
		sam->rnext[j++]=line[i++];
	}

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			temp[j]='\0';
			break;			
		}
		temp[j++]=line[i++];
	}
	sam->pnext=atoi(temp);

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			temp[j]='\0';
			break;			
		}
		temp[j++]=line[i++];
	}
	sam->tlen=atoi(temp);

	i++;
	j=0;
	sam->ncount=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->seq[j]='\0';
			break;			
		}
		if(line[i]=='N' || line[i]=='n')
		{
			sam->ncount++;
		}
		sam->seq[j++]=line[i++];
	}
	sam->readLength=j;

	i++;
	j=0;
	while(line[i])
	{
		if(line[i]=='\t' || line[i]==' ')
		{
			sam->qual[j]='\0';
			break;			
		}
		sam->qual[j++]=line[i++];
	}
	
	i++;
	
	char *temp2;

	temp2=strtok(&line[i],"\t\n ");
	
	while(temp2!=NULL)
	{
		if(temp2[0]=='M' && temp2[1]=='D')
		{
			strcpy(sam->md,temp2);
			break;
		}
		temp2=strtok(NULL,"\t\n ");
	}


	return sam;
}



void printHelp()
{

	cout<<"cgal v0.9.9-beta"<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"bowtie2convert - converts the map file in sam format outputted by Bowtie 2 into an internal format"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"bowtie2convert [options] <mapfile.sam> [maxFragmentLength]"<<endl; 
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<mapfile.sam>\t\t Map file outputted by Bowtie 2 in sam format"<<endl;
	cout<<endl;
	cout<<"Optional arguments:"<<endl;
	cout<<"[maxFragmentLength]\t Maximum insert size (fragment length). Default 5000"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<endl;
	exit(1);

}

int main(int argc, char *argv[])
{
	/*	input contig file name, read file name
		contig file - fasta format
		read file - fastq format		
	*/



	if(argc<2)
		printHelp();

	if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
		printHelp();


	if(argc==3)
		MAX_FRAGMENT_SIZE=atoi(argv[2]);
	else
		MAX_FRAGMENT_SIZE=5000;


	char *line= new char[MAX_REC_LEN];
	char *templine= new char[MAX_REC_LEN];
	
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	char * mapFileName=argv[1];


	mapFile=fopen(mapFileName, "r");

	if (mapFile == NULL) 
	{
		printf("Can't open map file\n");
		exit(1);
	}


	outFile=fopen("myout.sam","w");

	unFile=fopen("unmapped.txt","w");

	statFile=fopen("stat.txt","w");

	char *temp,nhstring[500];


	int it=0; 


	char preqname1[500];
	char preqname2[500];
	
	strcpy(preqname1,"*");
	strcpy(preqname2,"*");

	SAM *read1;
	SAM *read2;
		
	while(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
	{
		
		it++;

		
		if(line[0]=='@')
			continue;

		read1=getSAM(line);


		fgets(line, MAX_FILE_READ, mapFile);

		read2=getSAM(line);


		if(strcmp(read1->qname,preqname1)!=0 || strcmp(read2->qname,preqname2)!=0)
		{
			strcpy(preqname1,read1->qname);
			strcpy(preqname2,read2->qname);
			printVectors(outFile, unFile);
			reads1.push_back(read1);
			reads2.push_back(read2);
		}
		else
		{
			reads1.push_back(read1);
			reads2.push_back(read2);
		}	

	}

	

	fprintf(statFile,"%ld %ld %d %ld",totalCount, unCount, maxReadLength, MAX_FRAGMENT_SIZE*2);


	fclose(mapFile);
	fclose(outFile);
	fclose(unFile);
	fclose(statFile);

	


	return 0;
}
