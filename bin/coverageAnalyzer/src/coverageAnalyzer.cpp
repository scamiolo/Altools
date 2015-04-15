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
ofstream referenceCoverageFile;

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
int zeroCoverageWindows=0;
int	zeroCoverageWindowsReference =0;

float percentageCutoffLossAndGain;
float windowSize = 1;
float windowCoverage = 0;
float referenceWindowCoverage = 0;
float lossCutoff = 0;
float gainCutoff = 0;
float referenceAverageCoverage = 0;
float subjectAverageCoverage = 0;
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
	referenceCoverageFile.open("referenceCoverageFile");
	referenceCoverageFile << "Position	Chromosome	Coverage" << endl;


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
    cout.flush();

	system("Rscript DNAcopy.R");



	extractLossesAndGainsCutoff();

	cout << "Loss cutoff " << lossCutoff << endl;
	cout << "Gain cutoff " << gainCutoff << endl;

	generateLossesAndGainsOutput();

	system(("mv Rplots.pdf "+outputFolder).c_str());
	system(("mv cnv.txt " + outputFolder).c_str());
	system(("mv cnv_smoothed.txt " + outputFolder).c_str());
	system(("mv logRatioFile " + outputFolder).c_str());
	system(("mv referenceCoverageFile" + outputFolder).c_str());


	cout << "The reference average Coverage is " << referenceAverageCoverage << endl;
    cout.flush();
	cout << "The subject coverage average is " << subjectAverageCoverage << endl;
    cout.flush();
	getchar();
	getchar();

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

    cout.flush();
	//create the folder for the specific chromosome
	linuxCommand = "mkdir -p ";
	linuxCommand.append(outputFolder);
	linuxCommand.append(fileInInputFolder[chromosomeToAnalyze]);
	system(linuxCommand.c_str());

	//open zeroCoverage file
	

	//open zeroCoverage Gene Extractor file
	if(generateGEFiles==1)
		zeroCoverageGEFile.open( (outputFolder+fileInInputFolder[chromosomeToAnalyze]+"/"+
				fileInInputFolder[chromosomeToAnalyze]+"_zeroCoverage_GE").c_str() );


	//Load subject and reference genome coverage data in memory

	cout << "Loading data in memory......" << endl;
    cout.flush();
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

						if (fmod(analyzedLines ,3000000) == 0)
                        {
                            cout << analyzedLines << " analyzed in chromosome "
                            << fileInInputFolder[chromosomeToAnalyze] << endl;
                            cout.flush();
                        }
                        

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
    cout.flush();
	zeroCoverageFile << "Chromosome	Start	End	Length" << endl;

	zeroCoverageWindows=0;
	zeroCoverageWindowsReference =0;

	for (a=0;a<analyzedWindows;a++)
	{
		referenceCoverageFile << a*windowSize << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
				<< referenceCoverage[a] << endl;
		if(coverage[a]==0) zeroCoverageWindows++;
		else subjectAverageCoverage = subjectAverageCoverage + coverage[a];
		
		if(referenceCoverage[a]==0) zeroCoverageWindowsReference++;
		else referenceAverageCoverage = referenceAverageCoverage + referenceCoverage[a];
		
				
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
				<< deletionStart << "\t" << (int)(deletionStart+deletionArea) << "\t" << (int)deletionArea << endl;
		if (generateGEFiles == 1) zeroCoverageGEFile << "undefinedScaffold" << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
				<< deletionStart << "\t" << (int)(deletionStart+deletionArea+windowSize) << "\t" << (int)(deletionArea+windowSize) << endl;
		}
	}

	referenceAverageCoverage = referenceAverageCoverage / (analyzedWindows-1 - zeroCoverageWindowsReference);
	subjectAverageCoverage = subjectAverageCoverage / (analyzedWindows-1 - zeroCoverageWindows);
	




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
	ifstream referenceCoverageFileRead;

	ofstream losses;
	ofstream gains;
	ofstream losses_GE;
	ofstream gains_GE;

	ofstream losses_Smoothed;
	ofstream gains_Smoothed;
	ofstream losses_Smoothed_GE;
	ofstream gains_Smoothed_GE;
	
	ofstream testOutput;


	int a;
	int start;
	int end;
	int numberOfMarks;
	float segmentMean;
	float referencePosition;
	float startReference;
	float refCov;
	float segmentCoverageInReference;
	float numPositions;

	string header;
	char chromo[200];
	char previousChromo[200];
	char field[200];
	char filename_part[1000];
	char filename_Losses[1000];
	char filename_Gains[1000];
	char filename_Losses_GE[1000];
	char filename_Gains_GE[1000];
	char chromoReference[1000];
	char refCovChar[1000];
	char referencePositionChar[100];

	char filename_Smoothed_part[1000];
	char filename_Smoothed_Losses[1000];
	char filename_Smoothed_Gains[1000];
	char filename_Smoothed_Losses_GE[1000];
	char filename_Smoothed_Gains_GE[1000];


	testOutput.open("Test_output.txt");
	cnvFile.open("cnv.txt");
	referenceCoverageFile.close();
	referenceCoverageFileRead.open("referenceCoverageFile");


	//write losses and gain in each chromosome folder
	//read header from cnv file
	getline(cnvFile,header);
	getline(referenceCoverageFileRead,header);

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


		//Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
		referenceCoverageFileRead.seekg(0, ios::beg);
		getline(referenceCoverageFileRead,header);
		while(!referenceCoverageFileRead.eof())
		{
				referenceCoverageFileRead >> referencePosition;
				referenceCoverageFileRead >> chromoReference;
				referenceCoverageFileRead >> refCov;

				//look for chromosome in the reference coverage file
				//cout << chromo << "  "  << chromoReference << endl;
				//cout << strcmp(chromoReference,chromo) << endl;

				if(strcmp(chromoReference,chromo)==0)
				{
					//cout << "trovato chromo " << chromoReference << endl;
					//look for start position in reference coverage file
					while(referencePosition < start)
					{
						referenceCoverageFileRead >> referencePosition;
						referenceCoverageFileRead >> chromoReference;
						referenceCoverageFileRead >> refCov;
						if (referenceCoverageFileRead.eof())
						{
							segmentCoverageInReference = 0;
							break;
						}
					}
					numPositions = 0;
					segmentCoverageInReference = 0;
					while(referencePosition < end)
					{
						numPositions++;
						referenceCoverageFileRead >> referencePosition;
						referenceCoverageFileRead >> chromoReference;
						referenceCoverageFileRead >> refCov;
						segmentCoverageInReference = segmentCoverageInReference + refCov;
						if (referenceCoverageFileRead.eof())
						{
							segmentCoverageInReference = 0;
							break;
						}
					}
					segmentCoverageInReference = (segmentCoverageInReference/numPositions);
					break;
				}

		}

						testOutput << chromo << "  " << start << " SegMean: " << "SegMean: " << segmentMean << "first: " << segmentCoverageInReference / referenceAverageCoverage + 0.5 << " Second: " << subjectAverageCoverage/segmentCoverageInReference << endl;
						if((segmentMean)>= ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                gains << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								gains_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
							}


						if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                losses << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								losses_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end +windowSize - start) << endl;
							}


							//cout << "Ci arriva 1"  << endl;
		//Scan the rest of the file
		while(!cnvFile.eof())
		{
				for(a=0;a<2;a++)
				{
					cnvFile.get(field,100,' ');
					cnvFile.get();
				}
				if (cnvFile.eof()) break;
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


				//Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
				referenceCoverageFileRead.seekg(0, ios::beg);
				getline(referenceCoverageFileRead,header);
				while(!referenceCoverageFileRead.eof())
				{
						referenceCoverageFileRead >> referencePosition;
						referenceCoverageFileRead >> chromoReference;
						referenceCoverageFileRead >> refCov;

						//look for chromosome in the reference coverage file

						if(strcmp(chromoReference,chromo)==0)
						{
							cout << "trovato chromo " << chromoReference << endl;
							//look for start position in reference coverage file
							while(referencePosition < start)
							{
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
								if (referenceCoverageFileRead.eof())
								{
									segmentCoverageInReference = 0;
									break;
								}
							}
							numPositions = 0;
							segmentCoverageInReference = 0;
							while(referencePosition < end)
							{
								numPositions++;
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
								segmentCoverageInReference = segmentCoverageInReference + refCov;
								//cout << "refCov " << refCov << endl;
								if (referenceCoverageFileRead.eof())
								{
									segmentCoverageInReference = 0;
									break;
								}
							}
							segmentCoverageInReference = segmentCoverageInReference/numPositions;
							break;
						}

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

						testOutput << chromo << "  " << start << " SegMean: " << segmentMean << " first: " << segmentCoverageInReference / referenceAverageCoverage + 0.5 << " Second: " << subjectAverageCoverage/segmentCoverageInReference << endl;
						if((segmentMean)>=  ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                gains << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end +windowSize- start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								gains_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end +windowSize) << "\t" << (int)(end +windowSize - start) << endl;
							}


						if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                losses << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								losses_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end + windowSize - start) << endl;
							}
					}
					cout << "Ci arriva 2"  << endl;
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

				//Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
				referenceCoverageFileRead.seekg(0, ios::beg);
				getline(referenceCoverageFileRead,header);
				while(!referenceCoverageFileRead.eof())
				{
						referenceCoverageFileRead >> referencePosition;
						referenceCoverageFileRead >> chromoReference;
						referenceCoverageFileRead >> refCov;

						//look for chromosome in the reference coverage file

						if(strcmp(chromoReference,chromo)==0)
						{
							cout << "trovato chromo " << chromoReference << endl;
							//look for start position in reference coverage file
							while(referencePosition < start)
							{
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
									if (referenceCoverageFileRead.eof())
									{
										segmentCoverageInReference = 0;
										break;
									}
							}
							numPositions = 0;
							segmentCoverageInReference = 0;
							while(referencePosition < end)
							{
								numPositions++;
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
								segmentCoverageInReference = segmentCoverageInReference + refCov;
								//cout << "refCov " << refCov << endl;
								if (referenceCoverageFileRead.eof())
								{
									segmentCoverageInReference = 0;
									break;
								}
							}
							segmentCoverageInReference = segmentCoverageInReference/numPositions;
							break;
						}

				}

				cout << "Ci arriva 3"  << endl;
		if((segmentMean)>=  ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
			{
                gains_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
				if (generateGEFiles==1)
				gains_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
			}


		if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
			{
                losses_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
				if (generateGEFiles==1)
				losses_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
			}


		//Scan the rest of the file
		while(!cnvSmoothedFile.eof())
		{
				for(a=0;a<2;a++)
				{
					cnvSmoothedFile.get(field,100,' ');
					cnvSmoothedFile.get();
				}
				if (cnvSmoothedFile.eof()) break;
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


				//Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
				referenceCoverageFileRead.seekg(0, ios::beg);
				getline(referenceCoverageFileRead,header);
				while(!referenceCoverageFileRead.eof())
				{
						referenceCoverageFileRead >> referencePosition;
						referenceCoverageFileRead >> chromoReference;
						referenceCoverageFileRead >> refCov;

						//look for chromosome in the reference coverage file

						if(strcmp(chromoReference,chromo)==0)
						{
							cout << "trovato chromo " << chromoReference << endl;
							//look for start position in reference coverage file
							while(referencePosition < start)
							{
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
								if (referenceCoverageFileRead.eof())
								{
									segmentCoverageInReference = 0;
									break;
								}
							}
							numPositions = 0;
							segmentCoverageInReference = 0;
							while(referencePosition < end)
							{
								numPositions++;
								referenceCoverageFileRead >> referencePosition;
								referenceCoverageFileRead >> chromoReference;
								referenceCoverageFileRead >> refCov;
								segmentCoverageInReference = segmentCoverageInReference + refCov;
								//cout << "refCov " << refCov << endl;
								if (referenceCoverageFileRead.eof())
								{
									segmentCoverageInReference = 0;
									break;
								}
							}
							segmentCoverageInReference = segmentCoverageInReference/numPositions;
							break;
						}

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


						if(round(segmentMean)>=  round( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                gains_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								gains_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
							}


						if(round(segmentMean)<=  round( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(subjectAverageCoverage/segmentCoverageInReference) ) )
							{
                                losses_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                                 (segmentCoverageInReference/ referenceAverageCoverage) << "\t" << (segmentCoverageInReference*segmentMean)/subjectAverageCoverage << endl;
								if (generateGEFiles==1)
								losses_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
							}
					}
					//cout << "Ci arriva 4"  << endl;
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
