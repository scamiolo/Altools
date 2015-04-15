/*
 * pileupAnalyzer.cpp
 *
 *  Created on: Jan 27, 2011
 *      Author: salvo
 */

#include "pileupAnalyzer.h"
#include <iostream>
#include <fstream>



pileupAnalyzer::pileupAnalyzer() {
	// TODO Auto-generated constructor stub

}

void pileupAnalyzer::setPileupFolderName(string folderName) //Set the pileup folder path
{
	pileupFolderPath = folderName;
	return;
}

double pileupAnalyzer::getAverageCoverageInRange(string chromosome, int initialPosition, int finalPosition) //Get average coverage within a certain range
{
	double averageCoverage = 0;
	int a;
	int temp;
	int position, coverage, snp, indel;
	string prova;
	ifstream inputFile;

	if (initialPosition > finalPosition)
	{
		temp = initialPosition;
		initialPosition = finalPosition;
		finalPosition = temp;
	}
	inputFile.open( (pileupFolderPath + chromosome).c_str() );

	if(inputFile)
	{
		for(a=0 ; a<initialPosition-1; a++) inputFile >>position >> coverage >> snp >> indel;

		for(a=0 ; a<(finalPosition-initialPosition+1); a++)
		{

			inputFile >>position >> coverage >> snp >> indel;
			averageCoverage += coverage;
		}
		averageCoverage = averageCoverage / (finalPosition - initialPosition + 1);

	}
	else averageCoverage = -1;

	return averageCoverage;
}

int pileupAnalyzer::loadFiles(void) //Load files in the pileup folder
{
			inputFolder = opendir (pileupFolderPath.c_str());
			if (inputFolder != NULL) {

			/* load all the files and directories within directory */

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
			/* could not open directory */
			perror ("");
			return 0;//EXIT_FAILURE;
			}
			return 1;
}


void pileupAnalyzer::performSlidingAnalysis(string outputFolder, int windowSize, int step, int cutOff, int verbose)
{
	int a,i;
	int numberOfBases;
	int position;
	int *coverage = new int[1000000000];
	int *snp = new int[1000000000];
	int *indel = new int[1000000000];
	int totalCoverage;
	int totalIndel;
	int totalSnp;
	int windowPosition;
	int windowNumber;


	float overallAverageCoverage;
	float overallAverageIndel;
	float overallAverageSnp;
	float averageCoverage;
	float averageSnp;
	float averageIndel;

	ifstream inputFile;
	ofstream outputSlidingCoverage;
	ofstream outputSlidingSnp;
	ofstream outputSlidingIndel;


	for(a=1;a<numberOfFiles;a++)
	{
		inputFile.open( (pileupFolderPath + fileInInputFolder[a]).c_str() );
		cout << "Load chromosome " << fileInInputFolder[a] << " data in memory...." << endl;

		//Load data into memory
		numberOfBases = 0;
		totalCoverage = 0;
		totalSnp = 0;
		totalIndel = 0;

		while(!inputFile.eof())
		{
			numberOfBases++;
			inputFile >> position >> coverage[numberOfBases] >> snp[numberOfBases] >> indel[numberOfBases];
			if (verbose==1) cout << position << "  " << coverage[numberOfBases]<< "  " << snp[numberOfBases]<< "  " << indel[numberOfBases] << endl;
			if (coverage[numberOfBases]<cutOff) coverage[numberOfBases]=0;

			totalCoverage += coverage[numberOfBases];
			totalSnp += snp[numberOfBases];
			totalIndel += indel[numberOfBases];
		}

		//open output files
		outputSlidingCoverage.open( (outputFolder + fileInInputFolder[a] + "_CoverageWindows").c_str() );
		outputSlidingIndel.open( (outputFolder + fileInInputFolder[a] + "_snpWindows").c_str() );
		outputSlidingSnp.open( (outputFolder + fileInInputFolder[a] + "_IndelWindows").c_str() );

		//perform calculations
		cout << "Performing window calculations for chromosome " << fileInInputFolder[a] << endl;
		overallAverageCoverage = totalCoverage / (float)position;
		overallAverageSnp = totalSnp / (float)position;
		overallAverageIndel = totalIndel / (float)position;

		windowPosition = 1;
		windowNumber = 1;

		while(windowPosition <(position-windowSize))
		{
			averageCoverage = 0;
			averageSnp=0;
			averageIndel=0;

			for(i=windowPosition;i<(windowPosition + windowSize+1);i++)
			{
				averageCoverage = averageCoverage + coverage[i];
				averageSnp= averageSnp + snp[i];
				averageIndel= averageIndel + indel[i];
			}

			averageCoverage = averageCoverage/(float)windowSize;
			averageSnp = averageSnp/(float)windowSize;
			averageIndel = averageIndel/(float)windowSize;

			outputSlidingCoverage << windowNumber << "\t" << averageCoverage << "\t" << overallAverageCoverage << endl;
			outputSlidingIndel << windowNumber << "\t" << averageIndel << "\t" << overallAverageIndel << endl;
			outputSlidingSnp << windowNumber << "\t" << averageSnp << "\t" << overallAverageSnp << endl;

			windowNumber++;
			windowPosition=windowPosition + step;
		}

		inputFile.close();
		outputSlidingCoverage.close();
		outputSlidingIndel.close();
		outputSlidingSnp.close();
	}
}



pileupAnalyzer::~pileupAnalyzer() {
	// TODO Auto-generated destructor stub
}
