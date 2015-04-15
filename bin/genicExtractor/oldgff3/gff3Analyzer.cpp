/*
 * gff3Analyzer.cpp
 *
 *  Created on: Sep 9, 2011
 *      Author: salvo
 */

#include "gff3Analyzer.h"
#include <iostream>
#include <fstream>


string gff3FileName;

readAnnotationLine line[1000];
readAnnotationLine lineSorted[1000];
geneInformation *annotatedGene = new geneInformation[200000]; //the program will work only for 200000 genes
exonInformation *annotatedExon = new exonInformation[1000000]; // and 1000000 exons
geneInformation gene;
int exonLines;
int analyzedLines;
int isoformIsPresent = 0;
int numberOfAnnotatedGenes = 0;
int numberOfAnnotatedExons = 0;

ifstream fileToAnalyze;



gff3Analyzer::gff3Analyzer() {
	// TODO Auto-generated constructor stub


}

void gff3Analyzer::initializeGeneStatistics(void)
{
	exonLines=1;
	analyzedLines=-1;
	gene.numberOf3UtrExons = -1;
	gene.numberOf5UtrExons = -1;
	gene.numberOfCdsExons = -1;
	if (isoformIsPresent!=1) gene.geneName = "";
	gene.mRNAName = "";
}

int gff3Analyzer::getNumberOfAnnotatedGenes(void) //retrieve number of annotated genes
{
	return numberOfAnnotatedGenes;
}

geneInformation *gff3Analyzer::getAnnotatedGenes(void)
{
	return annotatedGene;
}

void gff3Analyzer::setFileToAnalyze(string fileName) // set the file name to work with
{
	gff3FileName = fileName;
}

void gff3Analyzer::collectAnnotationData(void) //load in memory the genic region position information
{
	int a;

	fileToAnalyze.open((gff3FileName).c_str());
	ignoreCommentLines();
	cout << "\nCollecting gene coordinates from annotation file..........." << endl;


	while(!fileToAnalyze.eof())
	{
		numberOfAnnotatedGenes++;
		initializeGeneStatistics();
		annotatedGene[numberOfAnnotatedGenes] = getGeneInformation();
		if (annotatedGene[numberOfAnnotatedGenes].geneName == "ENDOFFILE") //end of file reached!!
		{
			numberOfAnnotatedGenes--;
			break;
		}
	}
	cout << "Aquisition complete!" << endl;
	cout << "Number of found annotated genes = " << numberOfAnnotatedGenes << endl << "\nSorting the annotated genes database......" << endl;
	sortAnnotatedGenes(0,numberOfAnnotatedGenes);
	cout << "Sorting complete!" << endl;
	fileToAnalyze.close();

	return;
}


void gff3Analyzer::sortAnnotatedGenes(int left, int right) //sort the retrieved annotated genes by their initial position - use Quick Sort Algorithm
{

    int i = left, j = right;

    geneInformation tmp;

    int pivot = annotatedGene[(left + right) / 2].geneStartPosition;



    /* partition */

    while (i <= j) {

          while (annotatedGene[i].geneStartPosition < pivot)

                i++;

          while (annotatedGene[j].geneStartPosition > pivot)

                j--;

          if (i <= j) {

                tmp = annotatedGene[i];

                annotatedGene[i] = annotatedGene[j];

                annotatedGene[j] = tmp;

                i++;

                j--;

          }

    };



    /* recursion */

    if (left < j)

    	sortAnnotatedGenes(left, j);

    if (i < right)

    	sortAnnotatedGenes(i, right);

}

void gff3Analyzer::getExonsInformation(void)
{
	int a;
	readAnnotationLine line;
	string geneID;

	fileToAnalyze.open((gff3FileName).c_str());


	while(!fileToAnalyze.eof())
	{
		if(fileToAnalyze.eof()) break;

		line = readLineFromAnnotationFile();
		if(line.chromosome == "") break;
		if (line.region=="mRNA") geneID = extractIDFromAttribute(line.attribute);

		if(line.region!="gene" && line.region!="mRNA")
		{
			numberOfAnnotatedExons++;
			annotatedExon[numberOfAnnotatedExons].mRNAName = geneID;
			annotatedExon[numberOfAnnotatedExons].strand = line.strand;
			annotatedExon[numberOfAnnotatedExons].region = line.region;
			annotatedExon[numberOfAnnotatedExons].chromosome = line.chromosome;
			annotatedExon[numberOfAnnotatedExons].exonStartPosition = line.initialPosition;
			annotatedExon[numberOfAnnotatedExons].exonEndPosition = line.finalPosition;
			if(line.initialPosition > line.finalPosition)
			{
				annotatedExon[numberOfAnnotatedExons].exonStartPosition = line.finalPosition;
				annotatedExon[numberOfAnnotatedExons].exonEndPosition = line.initialPosition;
			}
		}
	}
	fileToAnalyze.close();
}

exonInformation *gff3Analyzer::getAnnotatedExons(void)
{
	return annotatedExon;
}

int gff3Analyzer::getNumberOfAnnotatedExons(void)
{
	return numberOfAnnotatedExons;
}


string gff3Analyzer::extractIDFromAttribute(string attribute)
{
	int a;
	string ID;


	for((a=attribute.find("ID=")+3); a<attribute.length(); a++)
	{
		if (attribute.at(a) != ';') ID += attribute.at(a);
		else break;
	}


	return ID;
}


geneInformation gff3Analyzer::getGeneInformation(void) //extract blocks of data from the annotation file relative to each gene
{


	int positionInFile;
	int a = 0;


	//read line relative to the gene
	if(isoformIsPresent!=1)
	{
		line[0] = readLineFromAnnotationFile();
		if (line[0].chromosome == "")
		{
			gene.geneName = "ENDOFFILE";
			return gene;
		}
		gene.geneStartPosition = line[0].initialPosition;
		gene.geneEndPosition = line[0].finalPosition;
		for((a=line[0].attribute.find("ID=")+3); a<line[0].attribute.length(); a++)
		{
			if (line[0].attribute.at(a) != ';') gene.geneName += line[0].attribute.at(a);
			else break;
		}
	}


	//read line relative to mRNA
	line[1] = readLineFromAnnotationFile();
	gene.mRNAStartPosition = line[1].initialPosition;
	gene.mRNAEndPosition = line[1].finalPosition;
	gene.strand = line[1].strand;
	gene.chromosome = line[1].chromosome;
	for((a=line[1].attribute.find("ID=")+3); a<line[1].attribute.length(); a++)
	{
		if (line[1].attribute.at(a) != ';') gene.mRNAName += line[1].attribute.at(a);
		else break;
	}

	//read line relative to the exons
	while(line[exonLines].region != "gene")
	{
		exonLines++;
		positionInFile = fileToAnalyze.tellg();
		line[exonLines] = readLineFromAnnotationFile();
		if(fileToAnalyze.eof()) break;
		if(line[exonLines].region=="mRNA")
		{
			isoformIsPresent = 1;
			break;
		}
		else isoformIsPresent = 0;
	}
	fileToAnalyze.seekg(positionInFile,ios_base::beg);
	exonLines--;

	//Elaboration of genic coordinates
	sortLines();
	assignLinesToRegions();

	return gene;
}

void gff3Analyzer::assignLinesToRegions(void) //go through the annotation lines and assign them to the regions 3UTR, 5UTR or CDS
{
	int a;

	analyzedLines = 0;
	for(a=2;a<=exonLines;a++)
	{
		analyzedLines++;
		assignPossibleUndefinedUTR(a);
		if (line[a].region == "5UTR")
		{
			gene.numberOf5UtrExons++;
			gene.utr5ExonStartPosition[gene.numberOf5UtrExons] = line[a].initialPosition;
			gene.utr5ExonEndPosition[gene.numberOf5UtrExons] = line[a].finalPosition;
		}


		if (line[a].region == "3UTR")
		{
			gene.numberOf3UtrExons++;
			gene.utr3ExonStartPosition[gene.numberOf3UtrExons]=line[a].initialPosition;
			gene.utr3ExonEndPosition[gene.numberOf3UtrExons]=line[a].finalPosition;
		}



		if (line[a].region == "CDS")
		{
			gene.numberOfCdsExons++;
			gene.cdsExonStartPosition[gene.numberOfCdsExons]=line[a].initialPosition;
			gene.cdsExonEndPosition[gene.numberOfCdsExons]=line[a].finalPosition;
		}

	}
}

void gff3Analyzer::assignPossibleUndefinedUTR(int lineNumber)
{
	if (line[lineNumber].region == "UTR") //if UTR code is used and no information on 5UTR or 3UTR is available the program will try a guess
			{
				if (checkFollowingCDSLines(1+analyzedLines)==true) //check whether there is a CDS line after the UTR
				{
					if(line[lineNumber].strand=="+")
					{
						line[lineNumber].region = "5UTR";
					}
					else line[lineNumber].region = "3UTR";
				}
				else
				{
					if(line[lineNumber].strand=="+")
					{
						line[lineNumber].region = "3UTR";
					}
					else line[lineNumber].region = "5UTR";
				}
			}
}

bool gff3Analyzer::checkFollowingCDSLines(int linesAlreadyAnalyzed) //while analyzing the annotation lines relative to a gene, check whether the remaining lines feature a CDS
{
	int a;
	for(a=linesAlreadyAnalyzed;a<=exonLines;a++)
	{
		if(line[a].region=="CDS") return true;
	}
	return false;
}

void gff3Analyzer::sortLines(void) //sort the annotation lines in ascending order of the initial exon position - Use bubble sort algorithm
{
	int i, j, flag = 1;    // set flag to 1 to start first pass
	readAnnotationLine temp;             // holding variable
	     for(i = 2; (i <= exonLines) && flag; i++)
	     {
	          flag = 0;
	         for (j=2; j <= (exonLines -1); j++)
	         {
	               if (line[j+1].initialPosition < line[j].initialPosition)      // ascending order simply changes to <
	              {
	                    temp = line[j];             // swap elements
	                    line[j] = line[j+1];
	                    line[j+1] = temp;
	                    flag = 1;               // indicates that a swap occurred.
	               }
	          }
	     }
	     return;   //arrays are passed to functions by address; nothing is returned

}


void gff3Analyzer::ignoreCommentLines(void) //ignore lines with comment at the beginning of the annotation file
{
	int positionInFile;
	string commentLine = "generic string";

	positionInFile = fileToAnalyze.tellg();
	getline(fileToAnalyze,commentLine);

	while(commentLine.at(0)=='#' && commentLine.at(1)=='#')
		{
			positionInFile = fileToAnalyze.tellg();
			getline(fileToAnalyze,commentLine);
	}
	fileToAnalyze.seekg(positionInFile,ios_base::beg);
}


readAnnotationLine gff3Analyzer::readLineFromAnnotationFile(void) //read a line from the annotation file
{
	readAnnotationLine line;
	string string;
	int temp;


	getline(fileToAnalyze,line.chromosome,'\t');
	getline(fileToAnalyze,line.algorythm,'\t');
	getline(fileToAnalyze,line.region,'\t');
	getline(fileToAnalyze,string,'\t');
	line.initialPosition = atoi(string.c_str());
	getline(fileToAnalyze,string,'\t');
	line.finalPosition = atoi(string.c_str());

	//make sure that initial and final position are reported in the right order
	if(line.initialPosition > line.finalPosition)
	{
		temp=line.initialPosition;
		line.initialPosition = line.finalPosition;
		line.finalPosition = temp;
	}

	getline(fileToAnalyze,line.score,'\t');
	getline(fileToAnalyze,line.strand,'\t');
	getline(fileToAnalyze,line.phase,'\t');
	getline(fileToAnalyze,line.attribute,'\n');
	return line;
}

void gff3Analyzer::annotatedGeneStatistics(void)
{
	int a;
	cout << "Analizzo il gene " << annotatedGene[numberOfAnnotatedGenes].geneName << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].mRNAName << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].numberOfCdsExons+1 << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].numberOf3UtrExons+1 << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].numberOf5UtrExons+1 << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].geneStartPosition << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].geneEndPosition << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].mRNAStartPosition << endl;
	cout << annotatedGene[numberOfAnnotatedGenes].mRNAEndPosition << endl;
	cout << "Posizioni CDS : " << endl;
	for (a=0;a<=annotatedGene[numberOfAnnotatedGenes].numberOfCdsExons;a++) cout << "CDS exon " << a << " = " << annotatedGene[numberOfAnnotatedGenes].cdsExonStartPosition[a] << "  " << annotatedGene[numberOfAnnotatedGenes].cdsExonEndPosition[a] << endl;
	for (a=0;a<=annotatedGene[numberOfAnnotatedGenes].numberOf3UtrExons;a++) cout << "3utr exon " << a << " = " << annotatedGene[numberOfAnnotatedGenes].utr3ExonStartPosition[a] << "  " << annotatedGene[numberOfAnnotatedGenes].utr3ExonEndPosition[a] << endl;
	for (a=0;a<=annotatedGene[numberOfAnnotatedGenes].numberOf5UtrExons;a++) cout << "5utr exon " << a << " = " << annotatedGene[numberOfAnnotatedGenes].utr5ExonStartPosition[a] << "  " << annotatedGene[numberOfAnnotatedGenes].utr5ExonEndPosition[a] << endl;
	getchar();
}

gff3Analyzer::~gff3Analyzer() {
	// TODO Auto-generated destructor stub
}
