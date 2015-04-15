/*
 * blastWrapper.h
 *
 *  Created on: Oct 20, 2011
 *      Author: salvo
 */

#ifndef BLASTWRAPPER_H_
#define BLASTWRAPPER_H_

#include <string>

using namespace std;

struct blastData
{
	string query;
	int queryLength;
	int numberOfAlignment;
	string alignmentInSubject[100000];
	string querySequence[100000];
	string subjectSequence[100000];
	string consensusSequence[100000];
	string cigarCode[100000];
	int subjectLength[100000];
	int alignmentLength[100000];
	int queryStart[100000];
	int queryEnd[100000];
	int subjectStart[100000];
	int subjectEnd[100000];
	int numberOfAlignments;
	float score[100000];
	int bits[100000];
	float evalue[100000];
	int identity[100000];
	int positive[100000];
	int gaps[100000];
	int frame[100000];
	string queryStrand[100000];
	string subjectStrand[100000];
	string strand[100000];
	int numberOfMismatches[100000];
	int numberOfGaps[100000];
};

using namespace std;

class blastWrapper {
public:
	blastWrapper();
	virtual ~blastWrapper();
	void openFile(string);
	void getAlignment(void);
	void getAlignmentFormat6(void);
	void moveToDoubleNewLine(void);
	void readHeader(void);
	int returnEOF(void);
	void setBlastType(string);
	void convertToLastz(string filename);
	string generateCigarCode(string, string,string,string);

};

#endif /* BLASTWRAPPER_H_ */
