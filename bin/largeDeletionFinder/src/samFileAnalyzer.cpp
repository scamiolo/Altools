/*
 * samFileAnalyzer.cpp
 *
 *  Created on: Jan 25, 2011
 *      Author: salvo
 */

#include "samFileAnalyzer.h"
#include <stdlib.h>


string fileName;
ifstream inputFile;

samFileAnalyzer::samFileAnalyzer() {
	// TODO Auto-generated constructor stub

}

void samFileAnalyzer::setFileName(string filename)  //Set file name
{
	fileName = filename;
	return;
}



int samFileAnalyzer::openFile(void) //open a file format file
{
	inputFile.open(fileName.c_str());
	if(inputFile) return 1;
	else return 0;
}

int samFileAnalyzer::closeFile(void) //close a sam format file
{
	inputFile.close();
	if(inputFile) return 1;
	else return 0;
}

samFormat samFileAnalyzer::readSamFileLine(void) //read a line from the sam format file
{
	samFormat line;
	inputFile >> line.readName >> line.flag >> line.region >> line.position >> line.mapq
	>> line.cigar >> line.regionInMate >> line.matePosition >> line.inferredInsertSize
	>> line.sequence >> line.quality;
	getline(inputFile,line.optional,'\n');

	return line;
}

int samFileAnalyzer::check_eof(void) //check whether the end of file has been reached
{
	if (inputFile.eof()) return 1;
	else return 0;
}


string samFileAnalyzer::checkLongIndel(pairSamFormat linePair, int minShortInsert,
		int maxShortInsert, int minLongInsert, int maxLongInsert) //check whether a read pair
																//indicates a delation or a long insertion in the target genome
{
	string indel;
	int temp;
	string tempString;

	if (linePair.position1>linePair.position2)
	{
		temp = linePair.position1;
		linePair.position1 = linePair.position2;
		linePair.position2 = temp;
		tempString = linePair.sequence1;
		linePair.sequence1 = linePair.sequence2;
		linePair.sequence2 = temp;
	}
	if( ((linePair.flag1==97 && linePair.flag2 == 145) || (linePair.flag1==81 && linePair.flag2 == 161))
       && (linePair.region1 == linePair.region2) && ( (abs((int)linePair.position2 - (int)linePair.position1 + (int)linePair.sequence1.length()) >= minShortInsert)
                                                     && (abs((int)linePair.position2 - (int)linePair.position1 + (int)linePair.sequence1.length()) <= maxShortInsert) )
					&& linePair.cigar1!="*" && linePair.cigar2!="*" ) return "shorterInsert";

	if( ((linePair.flag1==97 && linePair.flag2 == 145) || (linePair.flag1==81 && linePair.flag2 == 161))
			&& (linePair.region1 == linePair.region2) && ( (abs((int)linePair.position2 - (int)linePair.position1+ (int)linePair.sequence1.length()) >= minLongInsert)
					&& (abs((int)linePair.position2 - (int)linePair.position1 + (int)linePair.sequence1.length()) <= maxLongInsert) )
					&& linePair.cigar1!="*" && linePair.cigar2!="*" ) return "longerInsert";
	return "none";
}



void samFileAnalyzer::sortIndelByRegion(indelFormat array[1000000], int arrayLength)
{
	   int i, j, flag = 1;    // set flag to 1 to start first pass
	   indelFormat temp;             // holding variable

	     for(i = 1; (i <= arrayLength) && flag; i++)
	     {
	          flag = 0;
	          for (j=0; j < (arrayLength -1); j++)
	         {
	               if (array[j+1].region > array[j].region)      // ascending order simply changes to <
	              {
	                    temp = array[j];             // swap elements
	                    array[j] = array[j+1];
	                    array[j+1]= temp;
	                    flag = 1;               // indicates that a swap occurred.
	               }
	          }
	     }
}

void samFileAnalyzer::sortIndelByInitialPosition(indelFormat array[1000000], int arrayLength)
{
	   int i, j, flag = 1;    // set flag to 1 to start first pass
	   indelFormat temp;             // holding variable

	     for(i = 1; (i <= arrayLength) && flag; i++)
	     {
	         flag = 0;
	         for (j=0; j < (arrayLength -1); j++)
	         {
	               if ( (array[j+1].position1 < array[j].position1) && (array[j+1].region == array[j].region))      // ascending order simply changes to <
	              {
	                    temp = array[j];             // swap elements
	                    array[j] = array[j+1];
	                    array[j+1]= temp;
	                    flag = 1;               // indicates that a swap occurred.
	               }
	          }
	     }
}



pairSamFormat samFileAnalyzer::readSamFileLinePair(void) //read a pair of reads in the sam format file
{
	pairSamFormat line;

	inputFile >> line.readName1 >> line.flag1 >> line.region1 >> line.position1 >> line.mapq1
	>> line.cigar1 >> line.regionInMate1 >> line.matePosition1 >> line.inferredInsertSize1
	>> line.sequence1 >> line.quality1;
	getline(inputFile,line.optional1,'\n');


	inputFile >> line.readName2 >> line.flag2 >> line.region2 >> line.position2 >> line.mapq2
	>> line.cigar2 >> line.regionInMate2 >> line.matePosition2 >> line.inferredInsertSize2
	>> line.sequence2 >> line.quality2;
	getline(inputFile,line.optional2,'\n');

	return line;
}

void samFileAnalyzer::readHeader(void) //read header from sam file
{
	int positionInFile;
	string headerLine;

	positionInFile = inputFile.tellg();
	getline(inputFile,headerLine);

	while(headerLine.at(0)=='@')
	{
		positionInFile = inputFile.tellg();
		getline(inputFile,headerLine);
	}

	inputFile.seekg(positionInFile, ios::beg);

	return;
}

samFileAnalyzer::~samFileAnalyzer() {
	// TODO Auto-generated destructor stub
}
