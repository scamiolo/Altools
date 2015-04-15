//============================================================================
// Name        : slidingAnalysis.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "stdlib.h"
#include <fstream>
#include <dirent.h>
#include "string.h"
#include "stdio.h"
#include <ctype.h>
#include <sys/types.h>




using namespace std;

string pileupFolderName;
string outputFolder;
string pileupFolderPath;

int a,i;
int numberOfBases;
int position;
int *coverage = new int[500000000];
int *snp = new int[500000000];
int *indel = new int[500000000];
int totalCoverage;
int totalIndel;
int totalSnp;
int windowPosition;
int windowNumber;
int coverageZeroInWindows;
int overallZeroCoverage;


float overallAverageCoverage;
float overallAverageIndel;
float overallAverageSnp;
float averageCoverage;
float averageSnp;
float averageIndel;

ifstream inputFile;
ofstream outputSliding;
ofstream totalCoverageFile;

int windowSize;
int windowStep;
int cutOff;
int verbose;
int numberOfFiles;
int	totalNumberOfBases = 0;
float totalGenomeCoverage = 0;
int	totalGenomeZeroCoverage = 0;

string fileInInputFolder[200];
char linuxCommand[1000];

DIR *inputFolder;
struct dirent *ent;

void setPileupFolderName(string);
int loadFiles(void);
void performSlidingAnalysis(string,int, int, int, int);



int main(int argc, char** argv) {


	pileupFolderName = argv[1];
	outputFolder = argv[2];
	windowSize = atoi(argv[3]);
	windowStep = atoi(argv[4]);
	cutOff = atoi(argv[5]);
	verbose = atoi(argv[6]);




	setPileupFolderName(pileupFolderName.c_str());

	loadFiles();

	//performSlidingAnalysis(outputFolder, windowSize, windowStep, cutOff, verbose);

	for(a=1;a<numberOfFiles+1;a++)
	{
		inputFile.open( (pileupFolderPath + fileInInputFolder[a]).c_str() );
		cout << "Load chromosome " << fileInInputFolder[a] << " data in memory...." << endl;

		//Load data into memory
		numberOfBases = 0;
		totalCoverage = 0;
		totalSnp = 0;
		totalIndel = 0;
		overallZeroCoverage = 0;

		while(!inputFile.eof())
		{
			numberOfBases++;
			
			inputFile >> position >> coverage[numberOfBases] >> snp[numberOfBases] >> indel[numberOfBases];
			if (verbose==1) cout << position << "  " << coverage[numberOfBases]<< "  " << snp[numberOfBases]<< "  " << indel[numberOfBases] << endl;
			if (coverage[numberOfBases]<cutOff) coverage[numberOfBases]=0;

			totalCoverage += coverage[numberOfBases];
			
			totalSnp += snp[numberOfBases];
			totalIndel += indel[numberOfBases];
			if(coverage[numberOfBases]==0)
			{
				overallZeroCoverage++;
				
			}
		}
		totalGenomeCoverage += totalCoverage;
		totalNumberOfBases +=  numberOfBases;
		totalGenomeZeroCoverage += overallZeroCoverage;
		
		cout << totalGenomeCoverage << "  "  << totalNumberOfBases << "  "  << totalGenomeZeroCoverage << endl;
		
		//open output files
		outputSliding.open("slidingFile" );
		outputSliding << (fileInInputFolder[a]+"_Window").c_str() << "\taverageCoverage\toverallAverageCoverage\taverageSNP\toverallAverageSNP\taverageIndel\toverallAverageIndel" << endl;

		//perform calculations
		cout << "Performing window calculations for chromosome " << fileInInputFolder[a] << endl;
		overallAverageCoverage = totalCoverage / ((float)numberOfBases-(float)overallZeroCoverage);
		overallAverageSnp = totalSnp /  ((float)numberOfBases-(float)overallZeroCoverage);
		overallAverageIndel = totalIndel /  ((float)numberOfBases-(float)overallZeroCoverage);

		windowPosition = 1;
		windowNumber = 1;

		while(windowPosition <(position-windowSize))
		{
			averageCoverage = 0;
			averageSnp=0;
			averageIndel=0;
			coverageZeroInWindows=0;

			for(i=windowPosition;i<(windowPosition + windowSize+1);i++)
			{
				averageCoverage = averageCoverage + coverage[i];
				averageSnp= averageSnp + snp[i];
				averageIndel= averageIndel + indel[i];
				if (coverage[i]==0) coverageZeroInWindows++;
			}

			if ((windowSize-coverageZeroInWindows)!=0)
			{
				averageCoverage = averageCoverage/((float)windowSize-(float)coverageZeroInWindows);
				averageSnp = averageSnp/((float)windowSize-(float)coverageZeroInWindows);
				averageIndel = averageIndel/((float)windowSize-(float)coverageZeroInWindows);
			}
			else
			{
				averageCoverage = 0;
				averageSnp = 0;
				averageIndel = 0;
			}
			//averageCoverage = averageCoverage/((float)windowSize-(float)coverageZeroInWindows);
			//averageSnp = averageSnp/((float)windowSize-(float)coverageZeroInWindows);
			//averageIndel = averageIndel/((float)windowSize-(float)coverageZeroInWindows);

			outputSliding << windowNumber << "\t" << averageCoverage << "\t" << overallAverageCoverage <<
					"\t" << averageIndel << "\t" << overallAverageIndel << "\t" << averageSnp << "\t" <<
					overallAverageSnp << endl;
			windowNumber++;
			windowPosition=windowPosition + windowStep;
		}
		


		inputFile.close();
		outputSliding.close();

		system("R --save <createPlot.R");
		//move files in the output folder


		strcpy(linuxCommand,"mv slidingFile ");
		strcat(linuxCommand,outputFolder.c_str());
		strcat(linuxCommand,(fileInInputFolder[a] + "_Windows\n").c_str());
		system(linuxCommand);



		strcpy(linuxCommand,"mv Rplots.pdf ");
		strcat(linuxCommand,outputFolder.c_str());
		strcat(linuxCommand,(fileInInputFolder[a] + "Plots.pdf").c_str());
		system(linuxCommand);
	}
	cout << totalGenomeCoverage << "  "  << totalNumberOfBases << "  "  << totalGenomeZeroCoverage << endl;
		
	totalCoverageFile.open( (outputFolder + "totalCoverage").c_str() );
	totalGenomeCoverage = totalGenomeCoverage / ((float)totalNumberOfBases-(float)totalGenomeZeroCoverage);
	totalCoverageFile << "Genome average coverage: " << totalGenomeCoverage << endl;
	totalCoverageFile.close();
	
	cout << "\n\n\nProcess successfully completed!\nThanks for using Altools\nPress any key to continue.......\n\n" << endl;
	getchar();
	return 0;
}


void setPileupFolderName(string folderName) //Set the pileup folder path
{
	pileupFolderPath = folderName;
	return;
}


int loadFiles(void) //Load files in the pileup folder
{
			inputFolder = opendir (pileupFolderPath.c_str());
			if (inputFolder != NULL) {

			// load all the files and directories within directory

				numberOfFiles=0;
				while ((ent = readdir (inputFolder)) != NULL)
				{

						numberOfFiles++;
						fileInInputFolder[numberOfFiles] = ent->d_name;
						//extension = fileInInputFolder[numberOfFiles].substr(fileInInputFolder[numberOfFiles].length()-4, fileInInputFolder[numberOfFiles].length());
						if (fileInInputFolder[numberOfFiles] == "." || fileInInputFolder[numberOfFiles] == "..")  numberOfFiles--;
				}
			closedir (inputFolder);
			} else {
			// could not open directory
			perror ("");
			return 0;//EXIT_FAILURE;
			}
			return 1;
}


void performSlidingAnalysis(string outputFolder, int windowSize, int step, int cutOff, int verbose)
{



}
