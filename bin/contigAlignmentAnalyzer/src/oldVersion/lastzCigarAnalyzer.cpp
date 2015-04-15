/*
 * lastzCigarAnalyzer.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: salvo
 */

#include "lastzCigarAnalyzer.h"
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "completePileupAnalyzer.h"


string fileName;
string shortIndelFileName;
ofstream shortIndelFile;
ifstream fileToAnalyze;
ofstream structuralVariationFile;
cigarFormat alignment[10000];
int numberOfAlignments;
int performPileUpSearch = 0;
int insertionGE=0;
int deletionGE=0;
int translocationGE=0;
int inversionGE=0;
int transPlusInverGE = 0;
int duplicationGE = 0;

ofstream insertionsGE;
ofstream deletionsGE;
ofstream translocationsGE;
ofstream inversionsGE;
ofstream transPlusInversGE;
ofstream duplicationsGE;


completePileupAnalyzer myPileup;

lastzCigarAnalyzer::lastzCigarAnalyzer() {
	// TODO Auto-generated constructor stub

}

void lastzCigarAnalyzer::setPileupParameters(int performPileupSearch, string folderName)
{
	myPileup.setFolder(folderName);
	performPileUpSearch = performPileupSearch;
	if(performPileUpSearch==1) myPileup.collectPolymorphysms();
}

void lastzCigarAnalyzer::setFile(string filename) //set the file name to work with
{
	fileName = filename;
	return;
}

int lastzCigarAnalyzer::getNumberOfAlignments(void)
{
	return numberOfAlignments;
}

void lastzCigarAnalyzer::openFile(void) //open the file to work with
{
	fileToAnalyze.open(fileName.c_str());
	if(fileToAnalyze) return;
	else
	{
		cout << "The selected input file does not exit\nThe program will now exit\n\n" << endl;
		exit(-1);
	}

}

void lastzCigarAnalyzer::writeUpShortIndel(string chromosome,int startPosition, int endPosition, string cigarCode)
{
	char flag;
	int a,num,foundIndel, foundSnp;
	string number;

	//the alingment is on the forward orientation of the chromosome
	if(endPosition>startPosition)
	{
		a=0;
		while(cigarCode.at(a)!='\n')
		{
			a++;
			flag = cigarCode.at(a);
			a=a+2;

			number="";
			while(cigarCode.at(a)!=' ' && cigarCode.at(a)!='\n')
			{
				number += cigarCode.at(a);
				
				if (isdigit(cigarCode.at(a))) a++;
				else break;
				if(cigarCode.length()<a) break;
			}
			num = atoi(number.c_str());

			if(performPileUpSearch==1)
			{
				if(flag=='D' || flag == 'I')
				{
					foundIndel = myPileup.searchIndelInPileup(chromosome,startPosition+num);
					if (foundIndel == 1)
					{
						if(flag == 'D') shortIndelFile << chromosome << "\t" <<startPosition+num << "\tdeletion\t" << "1" << endl;
						if(flag == 'I') shortIndelFile << chromosome << "\t" <<startPosition+num  << "\tinsertion\t" << "1" << endl;
					}
					else
					{
						if(flag == 'D') shortIndelFile << chromosome << "\t" <<startPosition+num  << "\tdeletion\t" << "0" << endl;
						if(flag == 'I') shortIndelFile << chromosome << "\t" <<startPosition+num  << "\tinsertion\t" << "0" << endl;
					}

				}
				if(flag=='X')
				{
					foundSnp = myPileup.searchSnpInPileup(chromosome,startPosition+num);
					if (foundSnp == 1) shortIndelFile << chromosome << "\t" <<startPosition+num  << "\tsnp\t" << "1" << endl;
					else shortIndelFile << chromosome << "\t" <<startPosition+num  << "\tsnp\t" << "0" << endl;

				}
			}
		}
	}
	else //the alingment is on the reverse orientation of the chromosome
	{
		a=0;
		while(cigarCode.at(a)!='\n')
		{
			a++;
			flag = cigarCode.at(a);
			a=a+2;

			number="";
			while(cigarCode.at(a)!=' ' && cigarCode.at(a)!='\n')
			{
				number += cigarCode.at(a);
				a++;
			}
			num = atoi(number.c_str());

			if(performPileUpSearch==1)
			{
				if(flag=='D' || flag == 'I')
				{
					foundIndel = myPileup.searchIndelInPileup(chromosome,endPosition+num-1);
					if (foundIndel == 1)
					{
						if(flag == 'D') shortIndelFile << chromosome << "\t" <<endPosition+num-1 << "\tdeletion\t" << "1" << endl;
						if(flag == 'I') shortIndelFile << chromosome << "\t" <<endPosition+num-1  << "\tinsertion\t" << "1" << endl;
					}
					else
					{
						if(flag == 'D') shortIndelFile << chromosome << "\t" <<endPosition+num-1  << "\tdeletion\t" << "0" << endl;
						if(flag == 'I') shortIndelFile << chromosome << "\t" <<endPosition+num-1  << "\tinsertion\t" << "0" << endl;
					}

				}
				if(flag=='X')
				{
					foundSnp = myPileup.searchSnpInPileup(chromosome,endPosition+num-1);
					if (foundSnp == 1) shortIndelFile << chromosome << "\t" <<endPosition+num-1  << "\tsnp\t" << "1" << endl;
					else shortIndelFile << chromosome << "\t" <<endPosition+num-1  << "\tsnp\t" << "0" << endl;

				}
			}
		}
	}


}

int lastzCigarAnalyzer::checkEndOfFile(void) //check whether the end of file has been reached
{
	if(fileToAnalyze.eof()) return 1;
	else return 0;
}

void lastzCigarAnalyzer::elaborateAlignments(cigarFormat *alignments, int numOfAlignments, int minInsertionSize,
		int minDeletionSize, int minInversionSize, int minTranslocationSize,
		int minTranslocationPlusInversionSize, int minDuplicationSize) //elaborates the alignments for a scaffold and writeup structural variations in an output file
{
	int a,b;
	int duplicationSize;
	int deletionSize;
	int insertionSize;
	string structuralVariation;

	cout << "Alignments : " << numOfAlignments << endl;



	for (a=0;a<=numOfAlignments;a++)
	{
		if (performPileUpSearch==1) writeUpShortIndel(alignments[a].chromosome, alignments[a].chromosomeInitialPosition, alignments[a].chromosomeFinalPosition, (alignments[a].cigarCode+"\n"));

		/*
		for(b=a+1;b<numOfAlignments;b++)
		{
			structuralVariation = "None";
			duplicationSize = 0;
			insertionSize = 0;
			deletionSize = 0;

			//Analyzing possible duplication size
			if ( (alignments[b].scaffoldInitialPosition < alignments[a].scaffoldFinalPosition) &&
				 (alignments[b].scaffoldFinalPosition > alignments[a].scaffoldFinalPosition) )
				duplicationSize = alignments[a].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition;

			if( (alignments[b].scaffoldInitialPosition<=alignments[a].scaffoldFinalPosition) &&
				(alignments[b].scaffoldFinalPosition <=alignments[a].scaffoldFinalPosition)	)
				duplicationSize = alignments[b].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition;



			if (alignments[b].scaffoldInitialPosition >  alignments[a].scaffoldFinalPosition) //No overlap
			{
				if(alignments[a].chromosome ==  alignments[b].chromosome) //alignments on the same chromosome
				{
					if(alignments[a].chromosomeStrand ==  alignments[b].chromosomeStrand) //alignments in same orientation
					{
						if(alignments[b].chromosomeInitialPosition > alignments[a].chromosomeFinalPosition) //alignments are consecutive
						{
							//either a deletion or an insertion is occurring
							//Insertion
							if( ((alignments[b].scaffoldInitialPosition - alignments[a].scaffoldFinalPosition) -
								 (alignments[b].chromosomeInitialPosition - alignments[a].chromosomeFinalPosition))
								 >= minInsertionSize)
								{
									structuralVariation = "Insertion";
									insertionSize = (alignments[b].scaffoldInitialPosition - alignments[a].scaffoldFinalPosition) -
											 (alignments[b].chromosomeInitialPosition - alignments[a].chromosomeFinalPosition);
								}
							//Deletion
							if( ((alignments[b].chromosomeInitialPosition - alignments[a].chromosomeFinalPosition) -
							     (alignments[b].scaffoldInitialPosition - alignments[a].scaffoldFinalPosition))
								>= minDeletionSize )
								{
									structuralVariation = "Deletion";
									deletionSize = (alignments[b].chromosomeInitialPosition - alignments[a].chromosomeFinalPosition) -
										     (alignments[b].scaffoldInitialPosition - alignments[a].scaffoldFinalPosition);
								}
						}
						if(alignments[b].chromosomeFinalPosition < alignments[a].chromosomeInitialPosition) //Alignments are not consecutive a translocation occurred
						{
							if( (alignments[b].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition) >= minTranslocationSize)
								structuralVariation = "IntraTranslocation";
						}
					}
					else //algnments are in different orientation - an inversion occurred
					{
						if( (alignments[b].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition)>= minInversionSize)
							structuralVariation = "IntraInversion";
					}
				}
				else //alignments are in different chromosomes
				{
					if(alignments[a].chromosomeStrand ==  alignments[b].chromosomeStrand) //alignments in same orientation
					{
						if( (alignments[b].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition)>= minTranslocationSize)
							structuralVariation = "InterTranslocation";
					}
					else //alignments are in different orientations - a translocation + inversion occurred
					{
						if( (alignments[b].scaffoldFinalPosition - alignments[b].scaffoldInitialPosition)>= minTranslocationPlusInversionSize)
							structuralVariation = "InterTranslocationPlusInversion";
					}
				}
			}
			if(alignments[b].scaffoldInitialPosition < alignments[a].scaffoldFinalPosition) //A duplication + tranlocation/inversion occurred
			{
				if(alignments[a].chromosome ==  alignments[b].chromosome) //alignments on the same chromosome
				{
					if(alignments[a].chromosomeStrand ==  alignments[b].chromosomeStrand) //alignments in same orientation
					{
						//A duplication is followed by a translocation
						if(duplicationSize >= minDuplicationSize)
							structuralVariation = "IntraDuplicationPlusTranslocation";
					}
					else //alignment in opposite direction - the duplication event is followed by an inversion
					{
						if(duplicationSize >= minDuplicationSize)
							structuralVariation = "IntraDuplicationPlusInversion";
					}
				}
				else //alignments are in different chromosomes
				{
					if(alignments[a].chromosomeStrand ==  alignments[b].chromosomeStrand) //alignments in same orientation
					{
						//A duplication is followed by a translocation
						if(duplicationSize >= minDuplicationSize)
							structuralVariation = "InterDuplicationPlusTranslocation";
					}
					else //alignment in opposite direction - the duplication event is followed by an inversion
					{
						if(duplicationSize >= minDuplicationSize)
							structuralVariation = "InterDuplicationPlusInversion";
					}
				}
			}

			//write-up the result in an outputfile
			if(structuralVariation != "None")
			{
				structuralVariationFile << alignments[a].scaffold << "\t" << alignments[a].scaffoldInitialPosition << "\t" <<
					alignments[a].scaffoldFinalPosition << "\t" << alignments[a].chromosome << "\t" <<
					alignments[a].chromosomeInitialPosition << "\t" << alignments[a].chromosomeFinalPosition <<  "\t" <<
					alignments[a].chromosomeStrand << "\t" <<
					alignments[b].scaffoldInitialPosition << "\t" << alignments[b].scaffoldFinalPosition << "\t" <<
					alignments[b].chromosome << "\t" << alignments[b].chromosomeInitialPosition << "\t" <<
					alignments[b].chromosomeFinalPosition << "\t" << alignments[b].chromosomeStrand <<
					"\t" << insertionSize << "\t" << deletionSize << "\t" << duplicationSize << "\t" <<
					structuralVariation << endl;
				writeUpGEFiles(structuralVariation, alignment[a].scaffold,alignment[a].chromosome, alignment[a].chromosomeInitialPosition, alignment[a].chromosomeFinalPosition,
						alignment[b].chromosome, alignment[b].chromosomeInitialPosition, alignment[b].chromosomeFinalPosition);
			}


		}*/
	}
	writeUpShortIndel(alignments[numOfAlignments-1].chromosome, alignments[numOfAlignments-1].chromosomeInitialPosition, alignments[numOfAlignments-1].chromosomeFinalPosition,  (alignments[numOfAlignments-1].cigarCode+"\n"));


}

void lastzCigarAnalyzer::writeUpGEFiles(string structuralVariation, string scaffold, string chromosome1, int start1, int end1, string chromosome2, int start2, int end2)
{

	if (insertionGE==1 && structuralVariation=="Insertion") insertionsGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;
	if (deletionGE==1 && structuralVariation=="Deletion") deletionsGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;;
	if (translocationGE==1 && (structuralVariation=="IntraTranslocation" || structuralVariation=="InterTranslocation")) translocationsGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;;
	if (inversionGE==1 && structuralVariation=="IntraInversion") inversionsGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;;
	if (transPlusInverGE==1 && structuralVariation=="InterTranslocationPlusInversion") transPlusInversGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;;
	if (duplicationGE==1 &&  (structuralVariation=="IntraDuplicationPlusTranslocation" || structuralVariation=="IntraDuplicationPlusInversion" || structuralVariation=="InterDuplicationPlusTranslocation" || structuralVariation=="InterDuplicationPlusInversion") ) duplicationsGE << scaffold << "\t" << chromosome1 << "\t" << start1 << "\t" << end1 << "\t" << (end1-start1) << endl << scaffold << "\t" << chromosome2 << "\t" << start2 << "\t" << end2 << "\t" << (end2-start2) << endl;;
return;
}

void lastzCigarAnalyzer::openGEFiles(string outputFile, int insertion, int deletion, int translocation, int inversion, int transPlusInver, int duplication)
{
	if (insertion==1)  insertionsGE.open((outputFile+"_insertions_GE").c_str());
	if (deletion==1)  deletionsGE.open((outputFile+"_deletions_GE").c_str());
	if (translocation==1)  translocationsGE.open((outputFile+"_translocations_GE").c_str());
	if (inversion==1)  inversionsGE.open((outputFile+"_inversions_GE").c_str());
	if (transPlusInver==1)  transPlusInversGE.open((outputFile+"_transPlusInvers_GE").c_str());
	if (duplication==1)  duplicationsGE.open((outputFile+"_duplications_GE").c_str());
return;
}

void lastzCigarAnalyzer::setGEFlags(int insertion, int deletion, int translocation, int inversion, int transPlusInver, int duplication)
{

	insertionGE = insertion;
	deletionGE = deletion;
	translocationGE = translocation;
	inversionGE = inversion;
	transPlusInverGE = transPlusInver;
	duplicationGE = duplication;
}


void lastzCigarAnalyzer::sortAlignments(cigarFormat *alignments, int numOfAlignments) //sort alignments for a scaffold relative to the initial position - uses bubble sort
{
	int i, j, flag = 1;    // set flag to 1 to start first pass
		cigarFormat temp;             // holding variable
		     for(i = 1; (i < numOfAlignments) && flag; i++)
		     {
		         flag = 0;
		         for (j=0; j <= (numOfAlignments -2); j++)
		         {
		               if (alignments[j+1].scaffoldInitialPosition < alignments[j].scaffoldInitialPosition)      // ascending order simply changes to <
		              {
		                    temp = alignments[j];             // swap elements
		                    alignments[j] = alignments[j+1];
		                    alignments[j+1] = temp;
		                    flag = 1;               // indicates that a swap occurred.
		               }
		          }
		     }
		     return;   //arrays are passed to functions by address; nothing is returned
}
void lastzCigarAnalyzer::openStructuralVariationFile(string outputfile)
{
	structuralVariationFile.open(outputfile.c_str());
	return;
}

void lastzCigarAnalyzer::openShortIndelFile(string filename)
{
	shortIndelFile.open(filename.c_str());
	return;
}

cigarFormat *lastzCigarAnalyzer::getAlignments(void)
{
	return alignment;
}

void lastzCigarAnalyzer::setShortIndelFile(string filename)
{
	shortIndelFileName = filename;
}

void lastzCigarAnalyzer::getScaffoldAlignments(int bitCutoff)
{

	numberOfAlignments = 0;
	string previousScaffold;
	int positionInFile;
	int temp;

		//get first alignment for a certain scaffold
		begin:
			fileToAnalyze >> alignment[numberOfAlignments].format >> alignment[numberOfAlignments].scaffold >>
			alignment[numberOfAlignments].scaffoldInitialPosition >> alignment[numberOfAlignments].scaffoldFinalPosition >>
			alignment[numberOfAlignments].scaffoldStrand >> 	alignment[numberOfAlignments].chromosome >>
			alignment[numberOfAlignments].chromosomeInitialPosition >> 	alignment[numberOfAlignments].chromosomeFinalPosition >>
			alignment[numberOfAlignments].chromosomeStrand >> 	alignment[numberOfAlignments].score;
			getline(fileToAnalyze,alignment[numberOfAlignments].cigarCode,'\n');
			if (alignment[numberOfAlignments].scaffoldStrand == alignment[numberOfAlignments].chromosomeStrand) alignment[numberOfAlignments].chromosomeStrand = "+";
			else alignment[numberOfAlignments].chromosomeStrand = "-";



			if (alignment[numberOfAlignments].scaffoldFinalPosition < alignment[numberOfAlignments].scaffoldInitialPosition )
			{
				temp = alignment[numberOfAlignments].scaffoldInitialPosition;
				alignment[numberOfAlignments].scaffoldInitialPosition = alignment[numberOfAlignments].scaffoldFinalPosition;
				alignment[numberOfAlignments].scaffoldFinalPosition = temp;
			}
			previousScaffold = alignment[numberOfAlignments].scaffold;

			if ( (alignment[numberOfAlignments].scaffoldFinalPosition - alignment[numberOfAlignments].scaffoldInitialPosition) <= bitCutoff )
			goto begin;



		//get the following alignments relative to the same scaffold.
		while(alignment[numberOfAlignments].scaffold==previousScaffold  && !fileToAnalyze.eof())
		{
			numberOfAlignments++;

			positionInFile = fileToAnalyze.tellg();
			fileToAnalyze >> alignment[numberOfAlignments].format >> alignment[numberOfAlignments].scaffold >>
			alignment[numberOfAlignments].scaffoldInitialPosition >> alignment[numberOfAlignments].scaffoldFinalPosition >>
			alignment[numberOfAlignments].scaffoldStrand >> 	alignment[numberOfAlignments].chromosome >>
			alignment[numberOfAlignments].chromosomeInitialPosition >> 	alignment[numberOfAlignments].chromosomeFinalPosition >>
			alignment[numberOfAlignments].chromosomeStrand >> 	alignment[numberOfAlignments].score;
			getline(fileToAnalyze,alignment[numberOfAlignments].cigarCode,'\n');
			if (alignment[numberOfAlignments].scaffoldStrand == alignment[numberOfAlignments].chromosomeStrand) alignment[numberOfAlignments].chromosomeStrand = "+";
			else alignment[numberOfAlignments].chromosomeStrand = "-";

			if (alignment[numberOfAlignments].scaffoldFinalPosition < alignment[numberOfAlignments].scaffoldInitialPosition )
			{
				temp = alignment[numberOfAlignments].scaffoldInitialPosition;
				alignment[numberOfAlignments].scaffoldInitialPosition = alignment[numberOfAlignments].scaffoldFinalPosition;
				alignment[numberOfAlignments].scaffoldFinalPosition = temp;
			}

			if ( (alignment[numberOfAlignments].scaffoldFinalPosition - alignment[numberOfAlignments].scaffoldInitialPosition) <= bitCutoff )
				numberOfAlignments--;

		}
		if (!fileToAnalyze.eof()) fileToAnalyze.seekg(positionInFile,ios::beg);


}





void lastzCigarAnalyzer::closeFile(void) //close file
{
	fileToAnalyze.close();
	return;
}

lastzCigarAnalyzer::~lastzCigarAnalyzer() {
	// TODO Auto-generated destructor stub
}
