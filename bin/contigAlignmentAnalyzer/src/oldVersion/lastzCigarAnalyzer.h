/*
 * lastzCigarAnalyzer.h
 *
 *  Created on: Oct 3, 2011
 *      Author: salvo
 */

#ifndef LASTZCIGARANALYZER_H_
#define LASTZCIGARANALYZER_H_

#include <string>

using namespace std;

struct cigarFormat
{
	string format;
	string scaffold;
	int scaffoldInitialPosition;
	int scaffoldFinalPosition;
	string scaffoldStrand;
	string chromosome;
	int chromosomeInitialPosition;
	int chromosomeFinalPosition;
	string chromosomeStrand;
	int score;
	string  cigarCode;
};

class lastzCigarAnalyzer {
public:
	lastzCigarAnalyzer();
	virtual ~lastzCigarAnalyzer();
	cigarFormat readCigarLine(void);
	void setFile(string);
	void setShortIndelFile(string);
	void openShortIndelFile(string);
	void openFile(void);
	void closeFile(void);
	int checkEndOfFile(void);
	void getScaffoldAlignments(int);
	void openStructuralVariationFile(string);
	void elaborateAlignments(cigarFormat*, int,int,int,int,int,int,int);
	int getNumberOfAlignments(void);
	void sortAlignments(cigarFormat*,int);
	void writeUpShortIndel(string,int,int,string);
	void setPileupParameters(int,string);
	void openGEFiles(string, int, int, int, int, int, int);
	void setGEFlags(int,int,int,int,int,int);
	void writeUpGEFiles(string, string,string, int, int, string, int, int);
	cigarFormat *getAlignments(void);

};

#endif /* LASTZCIGARANALYZER_H_ */
