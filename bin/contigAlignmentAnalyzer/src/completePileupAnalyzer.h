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
    
    int (*snpPosition)[10000000] = new int[100][10000000];
    
    int (*indelPosition)[10000000] = new int[100][10000000];
    
    int numberOfSnp[1000];
    int numberOfIndel[1000];
    
    int numberOfChromosomes;
    int numberOfFiles;
    string fileInPileupFolder[1000];
    
};

#endif /* COMPLETEPILEUPANALYZER_H_ */
