/*
 * multiFastaAnalyzer.h
 *
 *  Created on: Sep 8, 2011
 *      Author: salvo
 */

#include <string>


#ifndef MULTIFASTAANALYZER_H_
#define MULTIFASTAANALYZER_H_

using namespace std;



class multiFastaAnalyzer {


public:
	multiFastaAnalyzer();
	void collectSequences(string);
	void setFileName(string);


	string retrieveSequence(string);
	string fileToAnalyze;

	string *sequence;
	string *locus;
	string *title;
	int *positionInFile;

	int numberOfSequences;


	virtual ~multiFastaAnalyzer();
};

#endif /* MULTIFASTAANALYZER_H_ */
