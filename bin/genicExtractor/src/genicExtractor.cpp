//============================================================================
// Name        : genicExtractor.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "gff3Analyzer.h"


using namespace std;

string inputFile;
string annotationFile;
string scaffold;
string chromosome;

int initialPosition;
int finalPosition;
int difference;

int numOfAnnotatedGenes;
int numOfAnnotatedExons;

gff3Analyzer myGff3;
geneInformation *annotatedGenes;
exonInformation *annotatedExons;
ifstream inFile;
ofstream outFile;

void checkGenePresence(int position1, int position2, string chr);
void checkExonPresence(int position1, int position2, string chr);


int main(int argc, char** argv) {

	inputFile = argv[1];
	annotationFile = argv[2];

	//retrieve annotation information
	cout << "Loading annotation file in memory......." << endl;
	myGff3.setFileToAnalyze(annotationFile);
	myGff3.collectAnnotationData();
	myGff3.getExonsInformation();
	numOfAnnotatedGenes = myGff3.getNumberOfAnnotatedGenes();
	numOfAnnotatedExons = myGff3.getNumberOfAnnotatedExons();
	annotatedGenes = myGff3.getAnnotatedGenes();
	annotatedExons = myGff3.getAnnotatedExons();

	//open intput/output file
	inFile.open(inputFile.c_str());
	outFile.open( (inputFile + "_genes").c_str() );

	//check whether inputfile entries contains genes or exons
	while(!inFile.eof())
	{
		inFile >> scaffold >> chromosome >> initialPosition >> finalPosition >> difference ;
		cout << "Analyzing range " << initialPosition << "-" << finalPosition << " in chromosome " << chromosome << endl,
		//check for genes presence
		checkGenePresence(initialPosition,finalPosition,chromosome);

		//check for exons presence
		checkExonPresence(initialPosition,finalPosition,chromosome);
	}



	cout << "Process successfully completed!\nThanks for using Altools.\nPress any key to continue....." << endl;

	getchar();


	return 0;
}

void checkGenePresence(int position1, int position2, string chr)
{
	int a;
	for(a=0;a<numOfAnnotatedGenes;a++)
	{
		if(chr==annotatedGenes[a].chromosome)
		{
			//gene in completly contained
			if( (position1 <= annotatedGenes[a].mRNAStartPosition) && (position2 >= annotatedGenes[a].mRNAEndPosition) )
				outFile << scaffold << "\t" << chromosome << "\t" << position1 << "\t" << position2 << "\t" <<
				annotatedGenes[a].mRNAName << "\tgene\t" << "1\t" << endl;

			//gene is partially contained
			if( ((position1 >= annotatedGenes[a].mRNAStartPosition) && (position1 <= annotatedGenes[a].mRNAEndPosition) &&
				 (position2 >= annotatedGenes[a].mRNAEndPosition)) ||
				 ((position2 >= annotatedGenes[a].mRNAStartPosition) && (position2 <= annotatedGenes[a].mRNAEndPosition) &&
				  (position1 <= annotatedGenes[a].mRNAStartPosition)) )
				outFile << scaffold << "\t" << chromosome << "\t" << position1 << "\t" << position2 << "\t" <<
				annotatedGenes[a].mRNAName << "\tgene\t" << "0\t" << endl;
		}
	}
}


void checkExonPresence(int position1, int position2, string chr)
{
	int a;
	for(a=0;a<numOfAnnotatedExons;a++)
	{
		if(chr==annotatedExons[a].chromosome)
		{
			//gene in completly contained
			if( (position1 <= annotatedExons[a].exonStartPosition) && (position2 >= annotatedExons[a].exonEndPosition) )
				outFile << scaffold << "\t" << chromosome << "\t" << position1 << "\t" << position2 << "\t" <<
				annotatedExons[a].mRNAName << "\t" <<annotatedExons[a].region <<  "\t1" << endl;

			//gene is partially contained
			if( ((position1 >= annotatedExons[a].exonStartPosition) && (position1 <= annotatedExons[a].exonEndPosition) &&
				 (position2 >= annotatedExons[a].exonEndPosition)) ||
				 ((position2 >= annotatedExons[a].exonStartPosition) && (position2 <= annotatedExons[a].exonEndPosition) &&
				  (position1 <= annotatedExons[a].exonStartPosition)) )
				outFile << scaffold << "\t" << chromosome << "\t" << position1 << "\t" << position2 << "\t" <<
				annotatedExons[a].mRNAName << "\t" <<annotatedExons[a].region <<  "\t0" << endl;
		}
	}
}
