//============================================================================
// Name        : coverageAnalyzer.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>


using namespace std;

ifstream inputFile;
ifstream referenceFile;
ofstream outpuFile;

ofstream highCoverageFile;
ofstream highCoverageGEFile;
ofstream zeroCoverageFile;
ofstream zeroCoverageGEFile;
ofstream logRatioFile;

string fileInInputFolder[100000];
string fileInReferenceFolder[10000];
string inputDirectory;
string referencePileupFolder;
string outputFolder;

char command[1000];

int numberOfFiles,a;
int chromosomeToAnalyze;
int minimumReadsNumberLowCoverage;
int minimumReadsNumberHighCoverage;
int minimumAreaToOutput;
int verboseOutput;
int generateGEFiles;
int analyzedLines = 0;
int analyzedWindows = 0;
int analyzedReferenceLines = 0;
int deletionStart;
int deletionArea = 0;
int zeroCoverageInWindow = 0;
int zeroCoverageInReferenceWindow = 0;

float percentageCutoffLossAndGain;
float windowSize = 1;
float windowCoverage = 0;
float referenceWindowCoverage = 0;
float lossCutoff = 0;
float gainCutoff = 0;
float segmentCoverage[1000000];



DIR *inputFolder;
DIR *referencePileupDirectory;
struct dirent *ent;

int loadInputFiles(void);
float *coverage = new float[500000000];
float *referenceCoverage = new float[500000000];
int *sortedCoverage = new int[500000000];
int *sortedReferenceCoverage = new int[500000000];
void reportZeroCoverageAreas(void);
void reportHighDepthAreas(int coverageCutoff);
void sortCoverageValues(int left, int right);
void extractLossesAndGainsCutoff(void);
void generateLossesAndGainsOutput(void);
int extractHighCoverageCutoff(int *coverageArray, int percentage, int linesAnalyzed);


int main(int argc, char** argv) {


	//get parameters from gui;
	inputDirectory = argv[1];
	minimumReadsNumberLowCoverage = atoi(argv[2]);
	percentageCutoffLossAndGain = atof(argv[3]);
	minimumAreaToOutput = atoi(argv[4]);
	verboseOutput = atoi(argv[5]);
	generateGEFiles = atoi(argv[6]);
	outputFolder = argv[7];
	referencePileupFolder = argv[8];
	windowSize = atof(argv[9]);


	//Collect files in the pileup folder

	loadInputFiles();
	logRatioFile.open("logRatioFile"); //needed for dnaCopy for copy number variation estimation later on
	logRatioFile << "Position	Chromosome	logRatio" << endl;


	/* Analyze chromosomes reported in the pileup folder */
	chromosomeToAnalyze = 0;
	while(chromosomeToAnalyze<numberOfFiles)
	{
		chromosomeToAnalyze++;
		if(fileInInputFolder[chromosomeToAnalyze]!="." && fileInInputFolder[chromosomeToAnalyze]!="..")
		{
		//open file

		inputFile.open((inputDirectory+fileInInputFolder[chromosomeToAnalyze]).c_str());

		referenceFile.open((referencePileupFolder+"/"+fileInInputFolder[chromosomeToAnalyze]).c_str());

		//analyze zero coverage areas and high depth areas
		reportZeroCoverageAreas(); //within this function the high depth areas will be searched too


		inputFile.close();
		}
	}



	logRatioFile.close();
	cout << "Extracting copy number variations......" << endl;

	system("R --save <DNAcopy.R");



	extractLossesAndGainsCutoff();

	cout << "Loss cutoff " << lossCutoff << endl;
	cout << "Gain cutoff " << gainCutoff << endl;

	generateLossesAndGainsOutput();

	system(("mv Rplots.pdf "+outputFolder).c_str());
	system(("mv cnv.txt " + outputFolder).c_str());
	system(("mv cnv_smoothed.txt " + outputFolder).c_str());

	return 0;
}



void extractLossesAndGainsCutoff(void)
{
	ifstream cnvFile;
	string header;
	char field[200];
	int a;
	int numberOfSegments;


	cnvFile.open("cnv.txt");


	// read header
	getline(cnvFile,header);

	//Read Values
	numberOfSegments = 0;
	while(!cnvFile.eof())
	{
		for(a=0;a<6;a++)
		{
			cnvFile.get(field,100,' ');
			cnvFile.get();
		}
		cnvFile.get(field,100,'\n');
		segmentCoverage[numberOfSegments] = atof(field);
		numberOfSegments++;
	}numberOfSegments--;


	sortCoverageValues(0,numberOfSegments);


	lossCutoff = segmentCoverage[(int)(percentageCutoffLossAndGain*numberOfSegments/100)];
	gainCutoff = segmentCoverage[numberOfSegments - (int)(percentageCutoffLossAndGain*numberOfSegments/100)];


}

void reportZeroCoverageAreas(void)
{
	string linuxCommand;
	int a;




	int position1, coverage1, snp1, indel1;
	int position2, coverage2, snp2, indel2;
	int zeroCoverageDetected;
	analyzedLines = 0;

	cout << "Please wait while extracting zero-coverage areas from chromosome " <<
			fileInInputFolder[chromosomeToAnalyze]   << endl;

	//create the folder for the specific chromosome
	linuxCommand = "mkdir -p ";
	linuxCommand.append(outputFolder);
	linuxCommand.append(fileInInputFolder[chromosomeToAnalyze]);
	system(linuxCommand.c_str());

	//open zeroCoverage file
	zeroCoverageFile.open( (outputFolder+fileInInputFolder[chromosomeToAnalyze]+"/"+
			fileInInputFolder[chromosomeToAnalyze]+"_zeroCoverage").c_str() );

	//open zeroCoverage Gene Extractor file
	if(generateGEFiles==1)
		zeroCoverageGEFile.open( (outputFolder+fileInInputFolder[chromosomeToAnalyze]+"/"+
				fileInInputFolder[chromosomeToAnalyze]+"_zeroCoverage_GE").c_str() );


	//Load subject and reference genome coverage data in memory

	cout << "Loading data in memory......" << endl;
	analyzedWindows = 0;
	analyzedLines = 0;
	while(!inputFile.eof())
	{
		windowCoverage = 0;
		referenceWindowCoverage = 0;

		if(!inputFile.eof() && !referenceFile.eof())
		{

			zeroCoverageInWindow = 0;
			zeroCoverageInReferenceWindow = 0;

			for(a=0;a<windowSize;a++)
					{
						analyzedLines++;

						if (fmod(analyzedLines ,3000000) == 0) cout << analyzedLines << " analyzed in chromosome "
										<< fileInInputFolder[chromosomeToAnalyze] << endl;

						if(!inputFile.eof()) inputFile >> position1 >> coverage1 >> snp1 >> indel1;
						else break;

						

						if(!referenceFile.eof()) referenceFile >> position2 >> coverage2 >> snp2 >> indel2;
						
						
						else break;


					
						windowCoverage = windowCoverage + coverage1;
						if (coverage1 == 0) zeroCoverageInWindow++;
						referenceWindowCoverage = referenceWindowCoverage + coverage2;
						if (coverage2 == 0) zeroCoverageInReferenceWindow++;
					}
					if(windowSize!=zeroCoverageInWindow) coverage[analyzedWindows] = windowCoverage / (windowSize - zeroCoverageInWindow);
					else coverage[analyzedWindows] = 0;

					if (windowSize!=zeroCoverageInReferenceWindow) referenceCoverage[analyzedWindows] = referenceWindowCoverage / (windowSize - zeroCoverageInReferenceWindow);
					else referenceCoverage[analyzedWindows] = 0;
		}
		else break;

		analyzedWindows++;
		

	}
	inputFile.close();
	referenceFile.close();


	//fillup zeroCoverage files
	cout << "Extracring zero coverage data......" << endl;
	zeroCoverageFile << "Chromosome	Start	End	Length" << endl;



	for (a=0;a<analyzedWindows;a++)
	{
		if(coverage[a]!=0 && referenceCoverage[a]!=0)
				logRatioFile << a*windowSize << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
				<< (coverage[a]/referenceCoverage[a]) << endl;

		deletionArea = 0;
		deletionStart = a*windowSize;
		if(coverage[a]==0 && referenceCoverage[a]!=0)
		{
			deletionArea = windowSize;
			while(coverage[a+1]==0 && referenceCoverage[a]!=0)
			{
				deletionArea = deletionArea + windowSize;
				a=a+1;
				if(a==analyzedWindows) break;
			}
		if(deletionArea >= minimumAreaToOutput) 	zeroCoverageFile << fileInInputFolder[chromosomeToAnalyze] << "\t"
				<< deletionStart << "\t" << (deletionStart+deletionArea) << "\t" << deletionArea << endl;
		if (generateGEFiles == 1) zeroCoverageGEFile << "undefinedScaffold" << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
				<< deletionStart << "\t" << (deletionStart+deletionArea) << "\t" << deletionArea << endl;
		}
	}





	 zeroCoverageFile.close();
	 highCoverageFile.close();
	 if( generateGEFiles==1 )
	 {
		 zeroCoverageGEFile.close();
		 highCoverageGEFile.close();
	 }
}


int loadInputFiles(void)
{

			inputFolder = opendir (inputDirectory.c_str());
			referencePileupDirectory = opendir (referencePileupFolder.c_str());
			if (inputFolder != NULL) {

			/* load all the files and directories within directory */

				numberOfFiles=0;

				while ((ent = readdir (inputFolder)) != NULL)
				{
					numberOfFiles++;
					fileInInputFolder[numberOfFiles] = ent->d_name;
					fileInReferenceFolder[numberOfFiles] = ent->d_name;
				}
			closedir (inputFolder);
			} else {
			/* could not open directory */
			perror ("");
			return EXIT_FAILURE;
			}
			return 1;
}



void sortCoverageValues(int left, int right)
{

	   int i = left, j = right;

	    float tmp;

	    float pivot = segmentCoverage[(left + right) / 2];



	    /* partition */

	    while (i <= j) {

	          while (segmentCoverage[i] < pivot)

	                i++;

	          while (segmentCoverage[j] > pivot)

	                j--;

	          if (i <= j) {

	                tmp = segmentCoverage[i];

	                segmentCoverage[i] = segmentCoverage[j];

	                segmentCoverage[j] = tmp;

	                i++;

	                j--;

	          }

	    };



	    /* recursion */

	    if (left < j)

	    	sortCoverageValues(left, j);

	    if (i < right)

	    	sortCoverageValues(i, right);
}








int extractHighCoverageCutoff(int *coverageArray, int percentage, int linesAnalyzed)
{
	double positionInArray = ((double)linesAnalyzed*95/100);

	return coverageArray[(int)positionInArray];
}



void generateLossesAndGainsOutput(void)
{

	ifstream cnvFile;
	ifstream cnvSmoothedFile;

	ofstream losses;
	ofstream gains;
	ofstream losses_GE;
	ofstream gains_GE;

	ofstream losses_Smoothed;
	ofstream gains_Smoothed;
	ofstream losses_Smoothed_GE;
	ofstream gains_Smoothed_GE;


	int a;
	int start;
	int end;
	int numberOfMarks;
	float segmentMean;

	string header;
	char chromo[200];
	char previousChromo[200];
	char field[200];
	char filename_part[1000];
	char filename_Losses[1000];
	char filename_Gains[1000];
	char filename_Losses_GE[1000];
	char filename_Gains_GE[1000];

	char filename_Smoothed_part[1000];
	char filename_Smoothed_Losses[1000];
	char filename_Smoothed_Gains[1000];
	char filename_Smoothed_Losses_GE[1000];
	char filename_Smoothed_Gains_GE[1000];


	cnvFile.open("cnv.txt");


	//write losses and gain in each chromosome folder
	//read header from cnv file
	getline(cnvFile,header);

	//get the first line

		for(a=0;a<2;a++)
		{
			cnvFile.get(field,100,' ');
			cnvFile.get();
		}

		cnvFile.get();
		cnvFile.get(chromo,100,'"');
		cnvFile.get();
		cnvFile.get();
		cnvFile.get(field,100,' ');
		cnvFile.get();
		start=atoi(field);
		cnvFile.get(field,100,' ');
		cnvFile.get();
		end=atoi(field);
		cnvFile.get(field,100,' ');
		cnvFile.get();
		numberOfMarks = atoi(field);
		cnvFile.get(field,100,'\n');
		cnvFile.get();
		segmentMean = atof(field);
		strcpy(previousChromo, chromo);
		strcpy(filename_part,outputFolder.c_str());
		strcat(filename_part,chromo);
		strcat(filename_part,"/");
		strcat(filename_part,chromo);

		strcpy(filename_Losses, filename_part);
		strcat(filename_Losses,"_Losses");
		strcpy(filename_Gains, filename_part);
		strcat(filename_Gains,"_Gains");

		if (generateGEFiles==1)
		{
			strcpy(filename_Losses_GE,filename_part);
			strcat(filename_Losses_GE,"_Losses_GE");
			strcpy(filename_Gains_GE,filename_part);
			strcat(filename_Gains_GE,"_Gains_GE");
		}


			losses.open(filename_Losses);
			gains.open(filename_Gains);
			if(generateGEFiles==1)
			{
				losses_GE.close();
				gains_GE.close();
				losses_GE.open(filename_Losses_GE);
				gains_GE.open(filename_Gains_GE);
			}


		if(segmentMean>=gainCutoff)
			{
                gains << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
				if (1==1) //Sostituire 1 con il GE_flag
				gains_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
			}


		if(segmentMean<=lossCutoff)
			{
                losses << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
				if (1==1) //Sostituire 1 con il GE_flag
				losses_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
			}



		//Scan the rest of the file
		while(!cnvFile.eof())
		{
				for(a=0;a<2;a++)
				{
					cnvFile.get(field,100,' ');
					cnvFile.get();
				}

				cnvFile.get();
				cnvFile.get(chromo,100,'"');
				cnvFile.get();
				cnvFile.get();
				cnvFile.get(field,100,' ');
				cnvFile.get();
				start=atoi(field);
				cnvFile.get(field,100,' ');
				cnvFile.get();
				end=atoi(field);
				cnvFile.get(field,100,' ');
				cnvFile.get();
				numberOfMarks = atoi(field);
				cnvFile.get(field,100,'\n');
				cnvFile.get();
				segmentMean = atof(field);

				strcpy(filename_part,outputFolder.c_str());
				strcat(filename_part,chromo);
				strcat(filename_part,"/");
				strcat(filename_part,chromo);

				strcpy(filename_Losses, filename_part);
				strcat(filename_Losses,"_Losses");
				strcpy(filename_Gains, filename_part);
				strcat(filename_Gains,"_Gains");



				if (generateGEFiles==1)
				{
					strcpy(filename_Losses_GE,filename_part);
					strcat(filename_Losses_GE,"_Losses_GE");
					strcpy(filename_Gains_GE,filename_part);
					strcat(filename_Gains_GE,"_Gains_GE");
				}


				if(!cnvFile.eof())
				{

						if(strcmp(previousChromo,chromo) != 0)
						{
							losses.close();
							gains.close();
							losses.open(filename_Losses);
							gains.open(filename_Gains);
							strcpy(previousChromo,chromo);
							if(generateGEFiles==1)
							{
								losses_GE.close();
								gains_GE.close();
								losses_GE.open(filename_Losses_GE);
								gains_GE.open(filename_Gains_GE);
							}
						}


						if(segmentMean>=gainCutoff)
							{
                                gains << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
								if (generateGEFiles==1)
								gains_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
							}


						if(segmentMean<=lossCutoff)
							{
                                losses << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
								if (generateGEFiles==1)
								losses_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
							}
						}
		}






	//write losses_Smoothed and gain from SMOOTHED file in each chromosome folder
	//read header from cnv file



	cnvSmoothedFile.open("cnv_smoothed.txt");
	getline(cnvSmoothedFile, header);


	//get the first line

		for(a=0;a<2;a++)
		{
			cnvSmoothedFile.get(field,100,' ');
			cnvSmoothedFile.get();
		}

		cnvSmoothedFile.get();
		cnvSmoothedFile.get(chromo,100,'"');

		cnvSmoothedFile.get();
		cnvSmoothedFile.get();
		cnvSmoothedFile.get(field,100,' ');
		cnvSmoothedFile.get();
		start=atoi(field);
		cnvSmoothedFile.get(field,100,' ');
		cnvSmoothedFile.get();
		end=atoi(field);
		cnvSmoothedFile.get(field,100,' ');
		cnvSmoothedFile.get();
		numberOfMarks = atoi(field);
		cnvSmoothedFile.get(field,100,'\n');
		cnvSmoothedFile.get();
		segmentMean = atof(field);
		strcpy(previousChromo, chromo);
		strcpy(filename_Smoothed_part,outputFolder.c_str());
		strcat(filename_Smoothed_part,chromo);
		strcat(filename_Smoothed_part,"/");
		strcat(filename_Smoothed_part,chromo);

		strcpy(filename_Smoothed_Losses, filename_Smoothed_part);
		strcat(filename_Smoothed_Losses,"_losses_Smoothed");
		strcpy(filename_Smoothed_Gains, filename_Smoothed_part);
		strcat(filename_Smoothed_Gains,"_gains_Smoothed");

		if (generateGEFiles==1)
		{
			strcpy(filename_Smoothed_Losses_GE,filename_Smoothed_part);
			strcat(filename_Smoothed_Losses_GE,"_losses_Smoothed_GE");
			strcpy(filename_Smoothed_Gains_GE,filename_Smoothed_part);
			strcat(filename_Smoothed_Gains_GE,"_gains_Smoothed_GE");
		}


			losses_Smoothed.open(filename_Smoothed_Losses);
			gains_Smoothed.open(filename_Smoothed_Gains);
			if(generateGEFiles==1)
			{
				losses_Smoothed_GE.close();
				gains_Smoothed_GE.close();
				losses_Smoothed_GE.open(filename_Smoothed_Losses_GE);
				gains_Smoothed_GE.open(filename_Smoothed_Gains_GE);
			}


		if(segmentMean>=gainCutoff)
			{
                gains_Smoothed << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
				if (generateGEFiles==1)
				gains_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
			}


		if(segmentMean<=lossCutoff)
			{
                losses_Smoothed << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
				if (generateGEFiles==1)
				losses_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
			}


		//Scan the rest of the file
		while(!cnvSmoothedFile.eof())
		{
				for(a=0;a<2;a++)
				{
					cnvSmoothedFile.get(field,100,' ');
					cnvSmoothedFile.get();
				}

				cnvSmoothedFile.get();
				cnvSmoothedFile.get(chromo,100,'"');
				cnvSmoothedFile.get();
				cnvSmoothedFile.get();
				cnvSmoothedFile.get(field,100,' ');
				cnvSmoothedFile.get();
				start=atoi(field);
				cnvSmoothedFile.get(field,100,' ');
				cnvSmoothedFile.get();
				end=atoi(field);
				cnvSmoothedFile.get(field,100,' ');
				cnvSmoothedFile.get();
				numberOfMarks = atoi(field);
				cnvSmoothedFile.get(field,100,'\n');
				cnvSmoothedFile.get();
				segmentMean = atof(field);

				strcpy(filename_Smoothed_part,outputFolder.c_str());
				strcat(filename_Smoothed_part,chromo);
				strcat(filename_Smoothed_part,"/");
				strcat(filename_Smoothed_part,chromo);

				strcpy(filename_Smoothed_Losses, filename_Smoothed_part);
				strcat(filename_Smoothed_Losses,"_losses_Smoothed");
				strcpy(filename_Smoothed_Gains, filename_Smoothed_part);
				strcat(filename_Smoothed_Gains,"_gains_Smoothed");



				if (generateGEFiles==1)
				{
					strcpy(filename_Smoothed_Losses_GE,filename_Smoothed_part);
					strcat(filename_Smoothed_Losses_GE,"_losses_Smoothed_GE");
					strcpy(filename_Smoothed_Gains_GE,filename_Smoothed_part);
					strcat(filename_Smoothed_Gains_GE,"_gains_Smoothed_GE");
				}

				if(!cnvSmoothedFile.eof())
				{

						if(strcmp(previousChromo,chromo) != 0)
						{
							losses_Smoothed.close();
							gains_Smoothed.close();
							losses_Smoothed.open(filename_Smoothed_Losses);
							gains_Smoothed.open(filename_Smoothed_Gains);
							strcpy(previousChromo,chromo);
							if(generateGEFiles==1)
							{
								losses_Smoothed_GE.close();
								gains_Smoothed_GE.close();
								losses_Smoothed_GE.open(filename_Smoothed_Losses_GE);
								gains_Smoothed_GE.open(filename_Smoothed_Gains_GE);
							}
						}


						if(segmentMean>=gainCutoff)
							{
                                gains_Smoothed << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
								if (generateGEFiles==1)
								gains_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
							}


						if(segmentMean<=lossCutoff)
							{
                                losses_Smoothed << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << "\t" << segmentMean << endl;
								if (generateGEFiles==1)
								losses_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << end << "\t" << (end - start) << endl;
							}
						}
		}


		cnvFile.close();
		cnvSmoothedFile.close();

		losses.close();
		gains.close();
		losses_GE.close();
		gains_GE.close();

		losses_Smoothed.close();
		gains_Smoothed.close();
		losses_Smoothed_GE.close();
		gains_Smoothed_GE.close();

}
