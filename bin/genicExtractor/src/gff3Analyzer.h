/*
 * gff3Analyzer.h
 *
 *  Created on: Sep 9, 2011
 *      Author: salvo
 */

#include <string>
#include <sstream>
#include <cstdlib>
#include <time.h>

#ifndef GFF3ANALYZER_H_
#define GFF3ANALYZER_H_

using namespace std;


struct readAnnotationLine
{
	string chromosome;
	string algorythm;
	string region;
	int initialPosition;
	int finalPosition;
	string score;
	string strand;
	string phase;
	string attribute;
};

struct exonInformation
{
	string mRNAName;
	int mRNAStartPosition;
	int mRNAEndPosition;
	int exonStartPosition;
	int exonEndPosition;
	string strand;
	string chromosome;
	string region;
};
struct geneInformation
{
	string geneName;
	string mRNAName;
	string strand;
	string chromosome;
	int cdsExonStartPosition[1000];
	int cdsExonEndPosition[1000];
	int numberOfCdsExons;
	int utr3ExonStartPosition[1000];
	int utr3ExonEndPosition[1000];
	int numberOf3UtrExons;
	int utr5ExonStartPosition[1000];
	int utr5ExonEndPosition[1000];
	int numberOf5UtrExons;
	int geneStartPosition;
	int geneEndPosition;
	int mRNAStartPosition;
	int mRNAEndPosition;
};




class gff3Analyzer {
public:
	gff3Analyzer();
	virtual ~gff3Analyzer();
	string extractIDFromAttribute(string);
	void setFileToAnalyze(string);
	void collectAnnotationData(void);
	void ignoreCommentLines(void);
	readAnnotationLine readLineFromAnnotationFile(void);
	void sortLines(void);
	void assignLinesToRegions(void);
	bool checkFollowingCDSLines(int);
	void assignPossibleUndefinedUTR(int);
	void annotatedGeneStatistics(void);
	void initializeGeneStatistics(void);
	void sortAnnotatedGenes(int, int);
	int getNumberOfAnnotatedGenes(void);
	int getNumberOfAnnotatedExons(void);
	geneInformation getGeneInformation(void);
	geneInformation *getAnnotatedGenes(void);
	exonInformation *getAnnotatedExons(void);
	void getExonsInformation(void);



};


#endif /* GFF3ANALYZER_H_ */
