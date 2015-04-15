//============================================================================
// Name        : largeDeletionFinder.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include <dirent.h>
#include "samFileAnalyzer.h"
#include "pileupAnalyzer.h"
#include "multiFastaAnalyzer.h"


using namespace std;

string inputDirectory;
string chromosome;
string fileInInputFolder[100];
string indel;
string outputFileSuffix;
string pileupFolder;
string sequenceFile;
string sequenceExtract;
string command;
string field;

char pos1a_toString[1000];
char pos2a_toString[1000];
char pos1b_toString[1000];
char pos2b_toString[1000];

DIR *inputFolder;
struct dirent *ent;

int numberOfFiles;
int checkInputFile;
int performCoverageAnalysis;
int minimumShortInsert, maximumShortInsert;
int minimumLongInsert, maximumLongInsert;
int numberOfTotalInsertions = 0;
int numberOfTotalDeletions = 0;
int distanceToCollapseIndel = 0;
int numberOfCollapsedDeletions = 0;
int numberOfCollapsedInsertions = 0;
int readsConfirmingInsertion[1000000];
int readsConfirmingDeletion[1000000];
int generateGEFile;
int matePairAnalyzed;
int minimumNumberOfReadsToConfirmIndel;
int maxCoverageInIndel;
int queryStart, queryEnd, refStart, refEnd;
int tempInt;
int value1, value2;
int startPoint, endPoint;
int a, temp;


double averageCoverage;
double score;

int loadInputFiles(void);
void sortDetectedIndel(void);
void collapseCloseInsertions(indelFormat array[1000000],  int arrayLength);
void collapseCloseDeletions(indelFormat array[1000000],  int arrayLength);
void sequencesAquisition(string fastaFile);

pairSamFormat matePair;
indelFormat totalInsertions[1000000];
indelFormat totalDeletions[1000000];
indelFormat collapsedInsertions[1000000];
indelFormat collapsedDeletions[1000000];
samFileAnalyzer mySamAnalyzer;


pileupAnalyzer myPileup;
multiFastaAnalyzer myFastaAnalyzer;

ofstream outputInsertionsFile;
ofstream outputDeletionsFile;
ofstream outputInsertionsFileGE;
ofstream outputDeletionsFileGE;
ofstream referenceExtract;
ofstream query;
ofstream blastElaborated_startDeletion;
ofstream blastElaborated_endDeletion;
ifstream unmappedReads;
ifstream outblast;
ifstream blastResult;


int main(int argc, char** argv) {

	inputDirectory = argv[1];
	minimumShortInsert = atoi(argv[2]);
	maximumShortInsert = atoi(argv[3]);
	minimumLongInsert = atoi(argv[4]);
	maximumLongInsert = atoi(argv[5]);
	distanceToCollapseIndel = atoi(argv[6]);
	performCoverageAnalysis = atoi(argv[7]);
	outputFileSuffix = argv[8];
	pileupFolder = argv[9];
	generateGEFile = atoi(argv[10]);
	minimumNumberOfReadsToConfirmIndel = atoi (argv[11]);
	maxCoverageInIndel = atoi(argv[12]);
    sequenceFile = argv[13];

	loadInputFiles();
    
    cout << "Collecting fasta file informations" << endl;
    sequencesAquisition(sequenceFile);
    
	matePairAnalyzed=0;
	//Check indel in loaded files
	for(a=1;a<=numberOfFiles;a++)
	{
			cout << (inputDirectory + fileInInputFolder[a]).c_str() << endl;
			mySamAnalyzer.setFileName( (inputDirectory + fileInInputFolder[a]).c_str() );
			checkInputFile = mySamAnalyzer.openFile();

			if(checkInputFile==1)
			{
				mySamAnalyzer.readHeader();
				while(mySamAnalyzer.check_eof()==0)
				{
					matePair = mySamAnalyzer.readSamFileLinePair();
					matePairAnalyzed++;
					if(fmod(matePairAnalyzed,100000)==0) cout << matePairAnalyzed << " pairs analyzed..." << endl;

					if(matePair.position1 > matePair.position2)
					{
						temp = matePair.position1;
						matePair.position1 = matePair.position2;
						matePair.position2 = temp;
					}

					indel = mySamAnalyzer.checkLongIndel(matePair, minimumShortInsert, maximumShortInsert,
							minimumLongInsert, maximumLongInsert);

					if (indel=="longerInsert" )
					{
						totalDeletions[numberOfTotalDeletions].region = matePair.region1;
						totalDeletions[numberOfTotalDeletions].position1 = matePair.position1 + matePair.sequence1.length();
						totalDeletions[numberOfTotalDeletions].position2 = matePair.position2;
						numberOfTotalDeletions++;
					}

					if (indel=="shorterInsert" )
					{
						totalInsertions[numberOfTotalInsertions].region = matePair.region1;
						totalInsertions[numberOfTotalInsertions].position1 = matePair.position1 + matePair.sequence1.length();
						totalInsertions[numberOfTotalInsertions].position2 = matePair.position2;
						numberOfTotalInsertions++;
					}
				}
			}
			mySamAnalyzer.closeFile();
		}

	//Elaborate insertions and deletions previously found
	sortDetectedIndel();
	cout << "Collapsing deletions...." << endl;
	collapseCloseDeletions(totalDeletions,numberOfTotalDeletions);
	cout << "Collapsing insertions...." << endl;
	collapseCloseInsertions(totalInsertions,numberOfTotalInsertions);

    
    
	//open up output files
	cout << "writing up output file....."  << endl;
	outputInsertionsFile.open( (inputDirectory+outputFileSuffix+"_Insertions").c_str() );
	outputDeletionsFile.open( (inputDirectory+outputFileSuffix+"_Deletions").c_str() );
	if(generateGEFile==1)
	{
		outputInsertionsFileGE.open( (inputDirectory+outputFileSuffix+"_Insertions_GE").c_str() );
		outputDeletionsFileGE.open( (inputDirectory+outputFileSuffix+"_Deletions_GE").c_str() );
	}


	//Perform eventual coverage analysis and output the results
	if (performCoverageAnalysis==1)
	{
		myPileup.setPileupFolderName(pileupFolder);

		for(a=1;a<numberOfCollapsedDeletions;a++)
		{
			if(readsConfirmingDeletion[a]>=minimumNumberOfReadsToConfirmIndel)
			{
				cout << "Extracting coverage for deletions in region " << collapsedDeletions[a].region
						<< " in range " << collapsedDeletions[a].position1 << " " << collapsedDeletions[a].position2 << endl;
				averageCoverage = myPileup.getAverageCoverageInRange(collapsedDeletions[a].region,
						collapsedDeletions[a].position1,collapsedDeletions[a].position2);
				if(averageCoverage<=maxCoverageInIndel)
				{
					outputDeletionsFile << collapsedDeletions[a].region << "\t" << collapsedDeletions[a].position1 << "\t"
							<< collapsedDeletions[a].position2 << "\t" << readsConfirmingDeletion[a] << "\t" << averageCoverage << endl;
				if(generateGEFile==1) outputDeletionsFileGE << "UndefinedScaffold\t" << collapsedDeletions[a].region << "\t" << collapsedDeletions[a].position1 << "\t"
						<< collapsedDeletions[a].position2 << "\t" << (collapsedDeletions[a].position2-collapsedDeletions[a].position1) << endl;
				}
			}
		}

		for(a=1;a<numberOfCollapsedInsertions;a++)
		{
			if(readsConfirmingInsertion[a]>=minimumNumberOfReadsToConfirmIndel)
			{
				cout << "Extracting coverage for insertions in region " << collapsedInsertions[a].region
						<< " in range " << collapsedInsertions[a].position1 << " " << collapsedInsertions[a].position2 << endl;
				averageCoverage = myPileup.getAverageCoverageInRange(collapsedInsertions[a].region,
						collapsedInsertions[a].position1,collapsedInsertions[a].position2);
				if(averageCoverage<=maxCoverageInIndel)
				{
					outputInsertionsFile << collapsedInsertions[a].region << "\t" << collapsedInsertions[a].position1 << "\t"
							<< collapsedInsertions[a].position2 << "\t" << readsConfirmingDeletion[a] << "\t" << averageCoverage << endl;
				if(generateGEFile==1) outputInsertionsFileGE << "UndefinedScaffold\t" << collapsedInsertions[a].region << "\t" << collapsedInsertions[a].position1 << "\t"
						<< collapsedInsertions[a].position2 << "\t" << (collapsedInsertions[a].position2-collapsedInsertions[a].position1) << endl;
				}
			}
		}
	}

	else
	{
		//only output indels
		cout << "Writing output......" << endl;
		for(a=1;a<numberOfCollapsedDeletions;a++)
		{
			if(readsConfirmingDeletion[a]>=minimumNumberOfReadsToConfirmIndel)
			{
				outputDeletionsFile << collapsedDeletions[a].region << "\t" << collapsedDeletions[a].position1 << "\t"
						<< collapsedDeletions[a].position2 << "\t" << readsConfirmingDeletion[a] << endl;
				if(generateGEFile==1) outputDeletionsFileGE << "UndefinedScaffold\t" << collapsedDeletions[a].region << "\t" << collapsedDeletions[a].position1 << "\t"
						<< collapsedDeletions[a].position2 << "\t" << (collapsedDeletions[a].position2-collapsedDeletions[a].position1) << endl;

			}
		}

		for(a=1;a<numberOfCollapsedInsertions;a++)
		{
			if(readsConfirmingInsertion[a]>=minimumNumberOfReadsToConfirmIndel)
			{
				outputInsertionsFile << collapsedInsertions[a].region << "\t" << collapsedInsertions[a].position1 << "\t"
						<< collapsedInsertions[a].position2 << "\t" << readsConfirmingInsertion[a] << endl;
				if(generateGEFile==1) outputInsertionsFileGE << "UndefinedScaffold\t" << collapsedInsertions[a].region << "\t" << collapsedInsertions[a].position1 << "\t"
						<< collapsedInsertions[a].position2 << "\t" << (collapsedInsertions[a].position2-collapsedInsertions[a].position1) << endl;

			}
		}
	}

	cout << "\n\nProcess terminated!\n Thanks for using Altools.\nPress any key to conitnue......... " << endl;
	getchar();
	return 0;
}

void sequencesAquisition(string fastaFile)
{
    myFastaAnalyzer.collectSequences(fastaFile);
}


void collapseCloseDeletions(indelFormat array[1000000], int arrayLength)
{
	int i;


	for(i=0;i<=arrayLength;i++)
	{
		numberOfCollapsedDeletions++;

		collapsedDeletions[numberOfCollapsedDeletions] = array[i];
		readsConfirmingDeletion[numberOfCollapsedDeletions] = 1;

		if(array[i].region == array[i+1].region)
		{
			while( (array[i+1].position1 - array[i].position1) <=  distanceToCollapseIndel && array[i].region == array[i+1].region)
			{
				if (array[i+1].position1 > array[i].position1) collapsedDeletions[numberOfCollapsedDeletions].position1 = array[i+1].position1;
				if (array[i+1].position2 < array[i].position2) collapsedDeletions[numberOfCollapsedDeletions].position2 = array[i+1].position2;
				i++;
				readsConfirmingDeletion[numberOfCollapsedDeletions]++;
			}
		}
	}
    

    
    //Code block to implement!
    
    
    cout << "Converting sam file into bam. Please wait....." << endl;
    command = "samtools view -bS " + inputDirectory + "*.sam  >output.bam";
    cout << command << endl;
    system(command.c_str());
    system("samtools sort output.bam output_sorted");
    system("samtools index output_sorted.bam");
    
    
    
    command = "cp " + sequenceFile + " ./refExtr124324";
    system(command.c_str());
    system("makeblastdb -in refExtr124324 -dbtype nucl");
    
    for (i=1;i<numberOfCollapsedDeletions;i++)
    {
        
        if(readsConfirmingDeletion[i]>=minimumNumberOfReadsToConfirmIndel)
        {
            cout << "Breakpoints are " << collapsedDeletions[i].position1 << " and " << collapsedDeletions[i].position2 << " " << (collapsedDeletions[i].position2 - collapsedDeletions[i].position1) <<endl;
            
            sprintf(pos1a_toString,"%d",collapsedDeletions[i].position1-1000);
            sprintf(pos1b_toString,"%d",collapsedDeletions[i].position1+1000);
            sprintf(pos2a_toString,"%d",collapsedDeletions[i].position2-1000);
            sprintf(pos2b_toString,"%d",collapsedDeletions[i].position2+1000);
            
            //cout << "Extracting unmapped reads with mapped mates in range " <<pos1_toString << "-" <<pos2_toString << endl;
            command = "";
            command = "samtools view output_sorted.bam " + collapsedDeletions[i].region + ":" + pos1a_toString + "-" + pos1b_toString + " > unmappedEnds.txt";
             system(command.c_str());
            command = "samtools view -f 4 -F 264 output_sorted.bam " + collapsedDeletions[i].region + ":" + pos1a_toString + "-" + pos1b_toString + " >> unmappedEnds.txt";
            system(command.c_str());
            command = "samtools view output_sorted.bam " + collapsedDeletions[i].region + ":" + pos2a_toString + "-" + pos2b_toString + " >> unmappedEnds.txt";
            system(command.c_str());
            command = "samtools view -f 4 -F 264 output_sorted.bam " + collapsedDeletions[i].region + ":" + pos2a_toString + "-" + pos2b_toString + " >> unmappedEnds.txt";
            system(command.c_str());
            
            
           unmappedReads.open("unmappedEnds.txt");
            query.open("query.fasta");
            while(!unmappedReads.eof())
            {
                for(a=0;a<7;a++) unmappedReads >> field;
                unmappedReads >> field;
                
                if (!unmappedReads.eof()) query << ">" << field << endl;
                unmappedReads >> field;
                unmappedReads >> field;
                
                if (!unmappedReads.eof()) query << field << endl;
                //unmappedReads >> field;
                getline(unmappedReads,field,'\n');
            }
            query.close();
            
            system("blastn -query query.fasta -db refExtr124324 -task blastn -outfmt 6 -out outblast.txt");
            //system("fasta36 query.fasta refExtr124324 -E 0.0001 -m 8 > fastaOutput");
            //system("sed 'd/#/' fastaOutput > outblast.txt");
            
            unmappedReads.close();
            //Elaborate blast result
            
            blastResult.open("outblast.txt");
            blastElaborated_startDeletion.open("elaboratedOutputBlast_startDeletion");
            blastElaborated_endDeletion.open("elaboratedOutputBlast_endDeletion");
            while(!blastResult.eof())
                  {
                      blastResult >> field >> chromosome >> score >> field >> field >>field >> queryStart >> queryEnd >> refStart >> refEnd >> field >> field;
                      if(refStart > refEnd)
                      {
                          tempInt = refStart;
                          refStart = refEnd;
                          refEnd = tempInt;
                      }
                      
                      //Change 70 with read length
                      if (chromosome==collapsedDeletions[i].region && (( abs(queryStart-queryEnd)<(int)(0.9*70))) && refStart <= collapsedDeletions[i].position2 && refEnd >= (collapsedDeletions[i].position1-70) && score>=97) //change with read length
                      {
                          //cout << refStart << " " << collapsedDeletions[i].position2 << " " << refEnd << " " << collapsedDeletions[i].position1 << endl;
                          
                          if( refEnd < (collapsedDeletions[i].position1+500) && refEnd>(collapsedDeletions[i].position1-500) )
                          blastElaborated_startDeletion << chromosome << "\t" << queryStart << "\t" << queryEnd << "\t" << refStart << "\t" << refEnd << endl;
                          if( refStart > (collapsedDeletions[i].position2-500) && refStart<(collapsedDeletions[i].position2+500) )
                          blastElaborated_endDeletion << chromosome << "\t" << queryStart << "\t" << queryEnd << "\t" << refStart << "\t" << refEnd << endl;
                      }

                  }
            blastResult.close();
            blastElaborated_startDeletion.close();
            blastElaborated_endDeletion.close();
            //get first breakpoint
            blastResult.open("elaboratedOutputBlast_startDeletion");
            startPoint = 0;
            while(!blastResult.eof())
            {
                blastResult >> field >> field >> field >> value1 >> value2;
                if(!blastResult.eof())
                {
                    if(value1 >= startPoint) startPoint = value1;
                    if(value2 >= startPoint) startPoint = value2;
                }
            }
            if(startPoint != 0 ) collapsedDeletions[i].position1 = startPoint+1;
            blastResult.close();
            
            //get second breakpoint
            blastResult.open("elaboratedOutputBlast_endDeletion");
            endPoint = 1000000000;
            while(!blastResult.eof())
            {
                blastResult >> field >> field >> field >> value1 >> value2;
                if(!blastResult.eof())
                {
                    if(value1 <= endPoint) endPoint = value1;
                    if(value2 <= endPoint) endPoint = value2;
                }
            }
            if(endPoint != 1000000000 ) collapsedDeletions[i].position2 = endPoint-1;
            blastResult.close();
            cout << " Refined Breakpoints are " << collapsedDeletions[i].position1 << " and " << collapsedDeletions[i].position2 << " " << (collapsedDeletions[i].position2 - collapsedDeletions[i].position1) <<endl;
            //getchar();
          }
    }
    
    
	return;
}


void collapseCloseInsertions(indelFormat array[1000000], int arrayLength)
{
	int i,temp;

	for(i=0;i<arrayLength;i++)
	{
		numberOfCollapsedInsertions++;

		collapsedInsertions[numberOfCollapsedInsertions] = array[i];
		readsConfirmingInsertion[numberOfCollapsedInsertions] = 1;

		if(array[i].region == array[i+1].region)
		{
			while( (array[i+1].position1 - array[i].position1) <=  distanceToCollapseIndel && array[i].region == array[i+1].region)
			{
				if (array[i+1].position1 > array[i].position1) collapsedInsertions[numberOfCollapsedInsertions].position1 = array[i+1].position1;
				if (array[i+1].position2 < array[i].position2) collapsedInsertions[numberOfCollapsedInsertions].position2 = array[i+1].position2;
				i++;
				readsConfirmingInsertion[numberOfCollapsedInsertions]++;
			}
		}
	}
	return;
}




void sortDetectedIndel(void) //sort insertion and deletion by the chromosome and by the first position
{
	cout << "Sorting insertions by region...." << endl;
	mySamAnalyzer.sortIndelByRegion(totalInsertions,numberOfTotalInsertions);
	cout << "Sorting insertions by position...." << endl;
	mySamAnalyzer.sortIndelByInitialPosition(totalInsertions,numberOfTotalInsertions);

	cout << "Sorting deletions by region...." << endl;
	mySamAnalyzer.sortIndelByRegion(totalDeletions,numberOfTotalDeletions);
	cout << "Sorting deletions by position...." << endl;
	mySamAnalyzer.sortIndelByInitialPosition(totalDeletions,numberOfTotalDeletions);
}

int loadInputFiles(void) //load input files from input folder.
{
	string extension;

			inputFolder = opendir (inputDirectory.c_str());
			if (inputFolder != NULL) {

			/* load all the files and directories within directory */

				numberOfFiles=0;
				while ((ent = readdir (inputFolder)) != NULL)
				{

						numberOfFiles++;
						fileInInputFolder[numberOfFiles] = ent->d_name;
						//extension = fileInInputFolder[numberOfFiles].substr(fileInInputFolder[numberOfFiles].length()-4, fileInInputFolder[numberOfFiles].length());
						if (fileInInputFolder[numberOfFiles] == "." || fileInInputFolder[numberOfFiles] == "..")  numberOfFiles--;

						//if ( fileInInputFolder[numberOfFiles].length()>4 &
						//		fileInInputFolder[numberOfFiles].substr(fileInInputFolder[numberOfFiles].length()-4, fileInInputFolder[numberOfFiles].length())!=".sam") numberOfFiles--;
				}
			closedir (inputFolder);
			} else {
			/* could not open directory */
			perror ("");
			return EXIT_FAILURE;
			}
			return 1;
}

