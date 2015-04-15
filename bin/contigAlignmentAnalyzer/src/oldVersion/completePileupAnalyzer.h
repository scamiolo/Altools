/*
 * completePileupAnalyzer.h
 *
 *  Created on: Oct 10, 2011
 *      Author: salvo
 */

#ifndef COMPLETEPILEUPANALYZER_H_
#define COMPLETEPILEUPANALYZER_H_
#include <string>

using namespace std;




class completePileupAnalyzer {
public:
	completePileupAnalyzer();
	virtual ~completePileupAnalyzer();
	void setFolder(string);
	void collectPolymorphysms(void);
	int loadPileupFiles(void);
	int searchSnpInPileup(string,int);
	int searchIndelInPileup(string,int);
};

#endif /* COMPLETEPILEUPANALYZER_H_ */
