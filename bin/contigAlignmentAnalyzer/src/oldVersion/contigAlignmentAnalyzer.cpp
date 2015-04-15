//============================================================================
// Name        : contigAlignmentAnalyzer.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "lastzCigarAnalyzer.h"
#include "blastWrapperFormat0.h"


using namespace std;

string inputFileName;
string outputFileName;
string pileupFolderName;
string alignmentFormat;
string freeSlot1;

int minimumBitCutOff;
int check_endOfFile;
int minimumInsertionSize;
int minimumDeletionSize;
int minimumTranslocationSize;
int minimumInversionSize;
int minimumTranslocationPlusInversionSize;
int minimumDuplicationSize;
int searchGoldenPolymorphysmSet;
int generateInsertionsGE;
int generateDeletionsGE;
int generateTranslocationGE;
int generateInversionGE;
int generateTranslocPlusInversion;
int generateDuplicationGE;


int a;


cigarFormat *cigarAlignments;
blastWrapper myBlast;


lastzCigarAnalyzer myCigarFile;

int main(int argc, char** argv) {

	inputFileName = argv[1];
	/* portion from beta version
	minimumBitCutOff = atoi(argv[2]);
	freeSlot1 = argv[3];
	minimumInsertionSize = atoi(argv[4]);
	minimumDeletionSize = atoi(argv[5]);
	minimumTranslocationSize = atoi(argv[6]);
	minimumInversionSize = atoi(argv[7]);
	minimumTranslocationPlusInversionSize = atoi(argv[8]);
	minimumDuplicationSize=atoi(argv[9]);
	searchGoldenPolymorphysmSet = atoi(argv[10]);*/
	pileupFolderName = argv[2]; //it is [11] in the beta version

	/*
	generateInsertionsGE= atoi(argv[12]);
	generateDeletionsGE = atoi(argv[13]);
	generateTranslocationGE = atoi(argv[14]);
	generateInversionGE = atoi(argv[15]);
	generateTranslocPlusInversion = atoi(argv[16]);
	generateDuplicationGE = atoi(argv[17]);
	alignmentFormat = argv[18];*/


		myBlast.convertToLastz(inputFileName);
		inputFileName+= ".cigar";

	/*if(alignmentFormat=="blastFormat0")
	{
		myBlast.convertToLastz(inputFileName);
		inputFileName+= ".cigar";
	}*/



	myCigarFile.setFile(inputFileName);
	myCigarFile.openFile();
	//myCigarFile.openStructuralVariationFile(inputFileName+"_StructuralVariation");
	myCigarFile.openShortIndelFile(inputFileName+"_goldenPolymorphisms");

	/*if(generateInsertionsGE==1 || generateDeletionsGE==1 || generateTranslocationGE==1 || generateInversionGE ==1 ||
			generateTranslocPlusInversion==1 || generateDuplicationGE == 1)
	{
		myCigarFile.openGEFiles(inputFileName, generateInsertionsGE,generateDeletionsGE,generateTranslocationGE,
						generateInversionGE, generateTranslocPlusInversion, generateDuplicationGE);
		myCigarFile.setGEFlags(generateInsertionsGE,generateDeletionsGE,generateTranslocationGE,
								generateInversionGE, generateTranslocPlusInversion, generateDuplicationGE);
	}*/



	//if(searchGoldenPolymorphysmSet==1)
	myCigarFile.setPileupParameters(1,pileupFolderName);



	while(check_endOfFile==0)
	{

		myCigarFile.getScaffoldAlignments(1); //(minimumBitCutOff)
		cigarAlignments = myCigarFile.getAlignments();
		cout << "Elaborating alignments for scaffold " << cigarAlignments[0].scaffold << endl;

		myCigarFile.sortAlignments(cigarAlignments, myCigarFile.getNumberOfAlignments());

		myCigarFile.elaborateAlignments(cigarAlignments,myCigarFile.getNumberOfAlignments(), 1,
				1,1, 1,
				1,1);
		//(cigarAlignments,myCigarFile.getNumberOfAlignments(), minimumInsertionSize,
		//minimumDeletionSize,minimumInversionSize, minimumTranslocationSize,
		//minimumTranslocationPlusInversionSize,minimumDuplicationSize)

		check_endOfFile = myCigarFile.checkEndOfFile();

	}


	system(("sort "+inputFileName+"_goldenPolymorphisms"+" -k1,1 -k2n,2n | uniq -f 1 > tempGolden ; mv tempGolden "+inputFileName+"_goldenPolymorphisms").c_str());

	cout << "Process successfully terminated!\nThanks for using Altolls.\nPress any key to continue....." << endl;
	getchar();
	return 0;
}
