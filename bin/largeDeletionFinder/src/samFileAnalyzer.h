/*
 * samFileAnalyzer.h
 *
 *  Created on: Jan 25, 2011
 *      Author: salvo
 */
#include <string>
#include <iostream>
#include <fstream>


#ifndef SAMFILEANALYZER_H_
#define SAMFILEANALYZER_H_

using namespace std;


struct indelFormat
{
	string region;
	int position1;
	int position2;
};




struct samFormat
{
	string readName;
	int flag;
	string region;
	int position;
	int mapq;
	string cigar;
	string regionInMate;
	int matePosition;
	int inferredInsertSize;
	string sequence;
	string quality;
	string optional;
};

struct pairSamFormat
{
	string readName1;
	int flag1;
	string region1;
	int position1;
	int mapq1;
	string cigar1;
	string regionInMate1;
	int matePosition1;
	int inferredInsertSize1;
	string sequence1;
	string quality1;
	string optional1;

	string readName2;
	int flag2;
	string region2;
	int position2;
	int mapq2;
	string cigar2;
	string regionInMate2;
	int matePosition2;
	int inferredInsertSize2;
	string sequence2;
	string quality2;
	string optional2;
};

class samFileAnalyzer {
public:


	samFileAnalyzer();
	virtual ~samFileAnalyzer();
	void setFileName(string);
	int openFile(void);
	int closeFile(void);
	samFormat readSamFileLine(void);
	pairSamFormat readSamFileLinePair(void);
	void readHeader(void);
	string checkLongIndel(pairSamFormat, int, int, int, int);
	int check_eof(void);
	void sortIndelByRegion(indelFormat[1000000], int);
	void sortIndelByInitialPosition(indelFormat[1000000], int);

};

#endif /* SAMFILEANALYZER_H_ */
