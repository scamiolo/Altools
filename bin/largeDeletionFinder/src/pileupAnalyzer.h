/*
 * pileupAnalyzer.h
 *
 *  Created on: Jan 27, 2011
 *      Author: salvo
 */
#include <string>
#include <dirent.h>

using namespace std;


#ifndef PILEUPANALYZER_H_
#define PILEUPANALYZER_H_

class pileupAnalyzer {
public:
	pileupAnalyzer();
	void setPileupFolderName(string);
	double getAverageCoverageInRange(string, int, int);
	virtual ~pileupAnalyzer();
	int loadFiles(void);
	void performSlidingAnalysis(string,int, int, int, int);

	string pileupFolderPath;

	int numberOfFiles;
	string fileInInputFolder[200];

	DIR *inputFolder;
	struct dirent *ent;
};

#endif /* PILEUPANALYZER_H_ */
