//============================================================================
// Name        : variantAnalyzer.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include "multiFastaAnalyzer.h"
#include "gff3Analyzer.h"


using namespace std;

ifstream inputFile;
ofstream outputFile;

int a, numberOfNucleotidesAroundPolymorphysm, verboseOutput, phase, elaboratedLines;
int minimumCoverage;

float snpPValue;
float indelPValue;

char *outputFileName, *inputFileName;
char nucleotideCode[4], forwardCodonTable[15][4], reverseCodonTable[15][4];
string codon[65], aminoacid[65], annotationFile, sequencesFolder, frameShift;
string retrievedSequence, nearNucleotidesString;
string codonFromReference, codonFromTarget;
string aminoacidFromReference;
string codon2aminoacid(string codon);

struct pileupLine
{
	string chromosome, referenceBase, consensusBase, calledBases, calledQualities, varFrequency, finalFields;
	int position, depth, reads1, reads2, strand1, strand2, quality1, quality2;
	double pvalue;

};

pileupLine readLine;

gff3Analyzer myGff3Analyzer;
multiFastaAnalyzer cdsAnalyzer, UTR3Analyzer, UTR5Analyzer;

void initializeVariables(void);
void sequencesAquisition(string folder);
void annotationAquisition(string file);
void getIndelData(void);
void getSNPData(void);
string retrieveNearNucleotides(string sequenceToAnalyze, int numberOfNucleotides, int polymorphysmPosition);
string getIndelFrameShift(string);
string getCodonFromReference(string sequenceFromReference, int position);
string getCodonFromTarget(string codonFromReference, char snp, int phase, string strand);
string codon2aminoacid(string codonToConvert);
string getAminoacidFromReference(string codonFromReference);
pileupLine readInputFileLine(void);


int main(int argc, char** argv) {


	//retrieve arguments from GUI
	inputFileName = argv[1];
	annotationFile = argv[2];
	sequencesFolder = argv[3];
	//numberOfNucleotidesAroundPolymorphysm = atoi(argv[4]);
	verboseOutput = atoi(argv[4]);
	minimumCoverage = atoi(argv[5]);
	snpPValue = atof(argv[6]);
	indelPValue = atof(argv[7]);

	numberOfNucleotidesAroundPolymorphysm = 0;


	initializeVariables();

	sequencesAquisition(sequencesFolder);
	cout << "\n\nAquisition sequences from file completed!" << endl;
	annotationAquisition(annotationFile);
	elaboratedLines=0;

	//Start scanning input file
	while(!inputFile.eof())
	{
		readLine = readInputFileLine();
		elaboratedLines++;
		if (verboseOutput == 1) cout << readLine.chromosome << "  " << readLine.position << endl;
		if(fmod(elaboratedLines,100000)==0) cout << elaboratedLines << " lines computed" << endl;
		if(inputFile.eof()) break;
		if ( (readLine.consensusBase.at(0) == '+' || readLine.consensusBase.at(0) == '-' || readLine.consensusBase.at(0) == '*')
				&& readLine.depth >= minimumCoverage && readLine.pvalue <= indelPValue) getIndelData();
		if ( readLine.referenceBase.at(0) != readLine.consensusBase.at(0) && (readLine.consensusBase.at(0) != '+' && readLine.consensusBase.at(0) != '-' && readLine.consensusBase.at(0) != '*')
				&& readLine.depth >= minimumCoverage && readLine.pvalue <= snpPValue ) getSNPData();
	}


	cout << "Process successfully terminated! Thanks for using Altools.\nPress any key to continue....." << endl;
	getchar();
	return 0;
}


void getSNPData(void)
{
	int i,j;
	int positionInSequence = 0;


	for(a=0;a<myGff3Analyzer.getNumberOfAnnotatedGenes();a++)
	{
		if(readLine.position<=myGff3Analyzer.getAnnotatedGenes()[a].geneEndPosition && readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].geneStartPosition && readLine.chromosome == myGff3Analyzer.getAnnotatedGenes()[a].chromosome)
		{
			//the SNP is within a gene! Let's check in which region it is!
			//Check whether is in 5UTR
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf5UtrExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf5UtrExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j] <= readLine.position)
						{
							if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j]+1);
							else  positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j]+1);
						}
					}
					retrievedSequence = UTR5Analyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);
					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;


					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);

					//output statistics
					if(retrievedSequence !="None")
					{
					outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "5UTR\tSNP\t" <<
							readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << "--" << "\t" <<
							readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
							((float)positionInSequence / (float)retrievedSequence.length())  << "\t--\t--\t--\t--\t--\t--\t--\t--\t"
							<< endl; //nearNucleotidesString << endl;
					}
				}
			}

			//Check whether is in 3UTR
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf3UtrExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf3UtrExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j] <= readLine.position)
						{
							if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j]+1);
							else  positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j]+1);
						}
					}

					retrievedSequence = UTR3Analyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);

					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;
					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);

					//output Statistics
					if(retrievedSequence !="None")
					{
					outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
					myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "3UTR\tSNP\t" <<
					readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
					myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << "--" << "\t" <<
					readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
					((float)positionInSequence / (float)retrievedSequence.length())  << "\t--\t--\t--\t--\t--\t--\t--\t--\t" <<
					endl; //		nearNucleotidesString << endl;
					}
				}
			}


			//Check whether is in CDS
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOfCdsExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOfCdsExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j] <= readLine.position)
							{
								if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j]+1);
								else positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j]+1);
							}
					}
					retrievedSequence = cdsAnalyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);




					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;

					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);

					if(retrievedSequence !="None")
					{
						codonFromReference = getCodonFromReference(retrievedSequence, positionInSequence);
						aminoacidFromReference = getAminoacidFromReference(codonFromReference);
						codonFromTarget = getCodonFromTarget(codonFromReference, readLine.consensusBase.at(0), phase, myGff3Analyzer.getAnnotatedGenes()[a].strand);

						//output Statistics
						outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
								myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "CDS\tSNP\t" <<
								readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
								myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << "--" << "\t" <<
								readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
								 ((float)positionInSequence / (float)retrievedSequence.length()) << "\t" << codonFromReference <<
								 "\t" << aminoacidFromReference << "\t" << codonFromTarget  << nearNucleotidesString << endl;

					}

				}
			}

		}
	}
}


string getCodonFromReference(string sequenceFromReference, int position)
{
	int numberOfCodons;
	string codon = "";


	numberOfCodons = position / 3;
	phase = position % 3;


	if(phase != 0)
	{
			codon.append(1,sequenceFromReference.at(numberOfCodons*3));
			codon.append(1,sequenceFromReference.at(numberOfCodons*3+1));
			codon.append(1,sequenceFromReference.at(numberOfCodons*3+2));
	}
	else
	{
		if(numberOfCodons>0)
		{
			codon.append(1,sequenceFromReference.at((numberOfCodons-1)*3));
			codon.append(1,sequenceFromReference.at((numberOfCodons-1)*3+1));
			codon.append(1,sequenceFromReference.at((numberOfCodons-1)*3+2));
		}
		else
		{
			codon.append(1,sequenceFromReference.at(numberOfCodons*3));
			codon.append(1,sequenceFromReference.at(numberOfCodons*3+1));
			codon.append(1,sequenceFromReference.at(numberOfCodons*3+2));
		}

	}



	//cout << "Codone di partenza " << codon << endl;

	return codon;
}

string getAminoacidFromReference(string codonFromReference)
{
	string Aminoacid= "";
	Aminoacid = codon2aminoacid(codonFromReference);

	//cout << "Amminoacido di partenza : " << Aminoacid << endl;

	return Aminoacid;
}

string getCodonFromTarget(string codonFromReference, char snp, int phase, string strand)
{
	string codonInTarget;
	string aminoacidInTarget;
	string codonString;
	string Aminoacid = "--";
	int i,j;



	if(strand.at(0)=='+')
	{
		for(i=0;i<=14;i++)
		{
			if(snp==forwardCodonTable[i][0])
			{
				for (j=1;j<=3;j++)
				{
					if(forwardCodonTable[i][j]!='$')
					{
						codonInTarget=codonFromReference;
						if (phase==0) codonInTarget.at(2) = forwardCodonTable[i][j];
						if (phase==2) codonInTarget.at(1) = forwardCodonTable[i][j];
						if (phase==1) codonInTarget.at(0) = forwardCodonTable[i][j];
						codonString.append(codonInTarget);
						codonString.append(1,'\t');
						aminoacidInTarget=codon2aminoacid(codonInTarget);
						codonString.append(aminoacidInTarget);
						codonString.append(1,'\t');
					}
					else codonString.append("--\t--\t");
				}
			}
		}
	}
	else
	{
		for(i=0;i<=14;i++)
		{
			if(snp==reverseCodonTable[i][0])
			{
				for (j=1;j<=3;j++)
				{
					if(forwardCodonTable[i][j]!='$')
					{
						codonInTarget=codonFromReference;

						if (phase==0) codonInTarget.at(2) = reverseCodonTable[i][j];
						if (phase==2) codonInTarget.at(1) = reverseCodonTable[i][j];
						if (phase==1) codonInTarget.at(0) = reverseCodonTable[i][j];

						codonString.append(codonInTarget);
						codonString.append(1,'\t');
						aminoacidInTarget=codon2aminoacid(codonInTarget);
						codonString.append(aminoacidInTarget);
						codonString.append(1,'\t');
					}
					else codonString.append("--\t--\t");
				}
			}
		}
	}

	return codonString;
}

string codon2aminoacid(string codonToConvert)
{
	string Aminoacid = "--";
	int a = 0;
	int noCodon = 1;

	for(a=1;a<=64;a++)
	{
		noCodon=1;
		if (codonToConvert == codon[a])
		{
			noCodon=0;
			break;
		}
	}
	if (noCodon==1) return "Unknown";

	Aminoacid = aminoacid[a];

	return Aminoacid;
}


void getIndelData(void) //check whether the indel is within a  genic region and report statistics
{
	int i,j;
	int positionInSequence = 0;


	for(a=0;a<myGff3Analyzer.getNumberOfAnnotatedGenes();a++)
	{
		if(readLine.position<=myGff3Analyzer.getAnnotatedGenes()[a].geneEndPosition && readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].geneStartPosition && readLine.chromosome == myGff3Analyzer.getAnnotatedGenes()[a].chromosome )
		{
			//the indel is within a gene! Let's check in which region it is!
			//Check whether is in 5UTR
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf5UtrExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf5UtrExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j] <= readLine.position)
						{
							if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j]+1);
							else  positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].utr5ExonStartPosition[j]+1);
						}
					}

					retrievedSequence = UTR5Analyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);
					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;
					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);
					//output statistics
					if(retrievedSequence !="None")
					{
					outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "5UTR\tINDEL\t" <<
							readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << "--" << "\t" <<
							readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
							((float)positionInSequence / (float)retrievedSequence.length())  << "\t--\t--\t--\t--\t--\t--\t--\t--\t" <<
							endl; //nearNucleotidesString << endl;
					}
				}
			}



			//Check whether is in 3UTR
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf3UtrExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOf3UtrExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j] <= readLine.position)
						{
							if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j]+1);
							else  positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].utr3ExonStartPosition[j]+1);
						}
					}

					retrievedSequence = UTR3Analyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);
					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;
					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);

					//output Statistics
					if(retrievedSequence !="None")
					{
					outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "3UTR\tINDEL\t" <<
							readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << "--" << "\t" <<
							readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
							((float)positionInSequence / (float)retrievedSequence.length())  << "\t--\t--\t--\t--\t--\t--\t--\t--\t" <<
							endl; //nearNucleotidesString << endl;
					}
				}
			}

			//Check whether is in CDS
			for(i=0;i<=myGff3Analyzer.getAnnotatedGenes()[a].numberOfCdsExons;i++)
			{
				if (readLine.position>= myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[i] && readLine.position<= myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[i])
				{

					//Determine position in sequence
					positionInSequence=0;
					for(j=0;j<=myGff3Analyzer.getAnnotatedGenes()[a].numberOfCdsExons;j++)
					{
						if (myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j] <= readLine.position)
							{
								if (readLine.position>=myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[j]) positionInSequence += (myGff3Analyzer.getAnnotatedGenes()[a].cdsExonEndPosition[j] - myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j]+1);
								else positionInSequence += (readLine.position - myGff3Analyzer.getAnnotatedGenes()[a].cdsExonStartPosition[j]+1);
							}
					}

					frameShift = getIndelFrameShift(readLine.consensusBase);
					retrievedSequence = cdsAnalyzer.retrieveSequence(myGff3Analyzer.getAnnotatedGenes()[a].mRNAName);
					if( myGff3Analyzer.getAnnotatedGenes()[a].strand == "-") positionInSequence = retrievedSequence.length() - positionInSequence + 1;
					//nearNucleotidesString =  retrieveNearNucleotides(retrievedSequence, numberOfNucleotidesAroundPolymorphysm, positionInSequence);

					//output Statistics
					if(retrievedSequence !="None")
					{
					outputFile << myGff3Analyzer.getAnnotatedGenes()[a].mRNAName << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].chromosome << "\t" << "CDS\tINDEL\t" <<
							readLine.referenceBase << "\t" << readLine.consensusBase << "\t" << readLine.depth << "\t" <<
							myGff3Analyzer.getAnnotatedGenes()[a].strand << "\t" << frameShift << "\t" <<
							readLine.position << "\t" << positionInSequence << "\t" << retrievedSequence.length() << "\t" <<
							 ((float)positionInSequence / (float)retrievedSequence.length()) << "\t--\t--\t--\t--\t--\t--\t--\t--\t" <<
							 endl; // nearNucleotidesString << endl;
					}
				}
			}

		}
	}
	return;
}


string retrieveNearNucleotides(string sequenceToAnalyze, int numberOfNucleotides, int polymorphysmPosition) // retrieve the nucleotides that surround the polymorphysm
{
	int a;
	string nucleotidesInRange = "";

	for(a=0;a<numberOfNucleotides;a++)
	{
		if( (polymorphysmPosition-numberOfNucleotides+a)>=0 )
		{
			nucleotidesInRange += sequenceToAnalyze.at(polymorphysmPosition-numberOfNucleotides+a);
			nucleotidesInRange.append("\t");
		}
		else nucleotidesInRange.append("--\t");
	}

	for(a=0;a<numberOfNucleotides;a++)
	{
		if ( (polymorphysmPosition+numberOfNucleotides-a)<sequenceToAnalyze.length())
		{
			nucleotidesInRange += sequenceToAnalyze.at(polymorphysmPosition+numberOfNucleotides-a);
			nucleotidesInRange.append("\t");
		}
		else nucleotidesInRange.append("--\t");
	}
	return nucleotidesInRange;
}


string getIndelFrameShift(string indelAnnotation) //extract the frame Shift information from the indel annotation stirng
{
	string firstAllele, secondAllele, frameShiftString;
	stringstream firstAlleleLength, secondAlleleLength;

	firstAllele = indelAnnotation.substr(0,indelAnnotation.find_first_of('/'));
	secondAllele = indelAnnotation.substr(indelAnnotation.find_first_of('/')+1, indelAnnotation.length());

	firstAlleleLength << (firstAllele.length()-1);
	secondAlleleLength << (secondAllele.length()-1);

	if(firstAllele.at(0)=='*' && secondAllele.at(0)!='*') frameShiftString = secondAllele.at(0) + secondAlleleLength.str();
	if(firstAllele.at(0)!='*' && secondAllele.at(0)=='*') frameShiftString = firstAllele.at(0) + firstAlleleLength.str();
	if(firstAllele.at(0)!='*' && secondAllele.at(0)!='*') frameShiftString = firstAllele.at(0) + firstAlleleLength.str() + "/" + secondAllele.at(0) + secondAlleleLength.str();
	return frameShiftString;
}





void annotationAquisition(string file)
{

	myGff3Analyzer.setFileToAnalyze(file);
	myGff3Analyzer.collectAnnotationData();
}

pileupLine readInputFileLine()
{
	pileupLine line;
	inputFile >> line.chromosome >> line.position >> line.referenceBase >> line.consensusBase >> line.reads1 >> line.reads2 >> line.varFrequency >> line.strand1 >> line.strand2 >>line.quality1 >> line.quality2 >> line.pvalue;
	getline(inputFile,line.finalFields, '\n');
	line.depth = line.reads1 + line.reads2;

	return line;
}


void sequencesAquisition(string folder)
{
		cdsAnalyzer.collectSequences(folder+"cds");
		UTR5Analyzer.collectSequences(folder+"5utr");
		UTR3Analyzer.collectSequences(folder+"3utr");
}


void initializeVariables(void)
{
	string str;
	//open input/output streams
	inputFile.open(inputFileName);
	getline(inputFile,str,'\n');
	outputFileName = inputFileName;
	strcat(outputFileName,"_statistics");
	outputFile.open(outputFileName);


	//write up output file header
	outputFile << "Locus	Chromosome	Region	Type	RefBase	ReadBase	Coverage	F/R	Frame_Shift	PositionInGenome	PositionInSequence	SequenceLength	RelativePosition	Ref_Codon	Ref_Aminoacid	1st_AlleleCodon	1stAllaleAminoAcid	2ndAlleleCodon	2ndAlleleAminoAcid	3rdAlleleCodon	3rdAlleleAminoAcid\t";
	for(a=0;a<numberOfNucleotidesAroundPolymorphysm;a++) outputFile << "base-" << (numberOfNucleotidesAroundPolymorphysm-a) << "\t";
	for(a=0;a<numberOfNucleotidesAroundPolymorphysm;a++) outputFile <<"base+" << (a+1) << "\t";
	outputFile << endl;


	//initialize codon/aminoacids table
	codon[1]= "TTT";    aminoacid[1]= "Phenylalanine";
	codon[2]= "TTC";    aminoacid[2]= "Phenylalanine";
	codon[3]= "TTA";    aminoacid[3]= "Leucine";
	codon[4]= "TTG";    aminoacid[4]= "Leucine";
	codon[5]= "CTT";    aminoacid[5]= "Leucine";
	codon[6]= "CTC";    aminoacid[6]= "Leucine";
	codon[7]= "CTA";    aminoacid[7]= "Leucine";
	codon[8]= "CTG";    aminoacid[8]= "Leucine";
	codon[9]= "ATT";    aminoacid[9]= "Isoleucine";
	codon[10]= "ATC";    aminoacid[10]= "Isoleucine";
	codon[11]= "ATA";    aminoacid[11]= "Isoleucine";
	codon[12]= "ATG";    aminoacid[12]= "Methionine";
	codon[13]= "GTT";    aminoacid[13]= "Valine";
	codon[14]= "GTC";    aminoacid[14]= "Valine";
	codon[15]= "GTA";    aminoacid[15]= "Valine";
	codon[16]= "GTG";    aminoacid[16]= "Valine";
	codon[17]= "TCT";    aminoacid[17]= "Serine";
	codon[18]= "TCC";    aminoacid[18]= "Serine";
	codon[19]= "TCA";    aminoacid[19]= "Serine";
	codon[20]= "TCG";    aminoacid[20]= "Serine";
	codon[21]= "CCT";    aminoacid[21]= "Proline";
	codon[22]= "CCC";    aminoacid[22]= "Proline";
	codon[23]= "CCA";    aminoacid[23]= "Proline";
	codon[24]= "CCG";    aminoacid[24]= "Proline";
	codon[25]= "ACT";    aminoacid[25]= "Threonine";
	codon[26]= "ACC";    aminoacid[26]= "Threonine";
	codon[27]= "ACA";    aminoacid[27]= "Threonine";
	codon[28]= "ACG";    aminoacid[28]= "Threonine";
	codon[29]= "GCT";    aminoacid[29]= "Alanine";
	codon[30]= "GCC";    aminoacid[30]= "Alanine";
	codon[31]= "GCA";    aminoacid[31]= "Alanine";
	codon[32]= "GCG";    aminoacid[32]= "Alanine";
	codon[33]= "TAT";    aminoacid[33]= "Tyrosine";
	codon[34]= "TAC";    aminoacid[34]= "Tyrosine";
	codon[35]= "TAA";    aminoacid[35]= "Stop";
	codon[36]= "TAG";    aminoacid[36]= "Stop";
	codon[37]= "CAT";    aminoacid[37]= "Histidine";
	codon[38]= "CAC";    aminoacid[38]= "Histidine";
	codon[39]= "CAA";    aminoacid[39]= "Glutamine";
	codon[40]= "CAG";    aminoacid[40]= "Glutamine";
	codon[41]= "AAT";    aminoacid[41]= "Asparagine";
	codon[42]= "AAC";    aminoacid[42]= "Asparagine";
	codon[43]= "AAA";    aminoacid[43]= "Lysine";
	codon[44]= "AAG";    aminoacid[44]= "Lysine";
	codon[45]= "GAT";    aminoacid[45]= "Aspartic_acid";
	codon[46]= "GAC";    aminoacid[46]= "Aspartic_acid";
	codon[47]= "GAA";    aminoacid[47]= "Glutamic_acid";
	codon[48]= "GAG";    aminoacid[48]= "Glutamic_acid";
	codon[49]= "TGT";    aminoacid[49]= "Cysteine";
	codon[50]= "TGC";    aminoacid[50]= "Cysteine";
	codon[51]= "TGA";    aminoacid[51]= "Stop";
	codon[52]= "TGG";    aminoacid[52]= "Tryptophan";
	codon[53]= "CGT";    aminoacid[53]= "Arginine";
	codon[54]= "CGC";    aminoacid[54]= "Arginine";
	codon[55]= "CGA";    aminoacid[55]= "Arginine";
	codon[56]= "CGG";    aminoacid[56]= "Arginine";
	codon[57]= "AGT";    aminoacid[57]= "Serine";
	codon[58]= "AGC";    aminoacid[58]= "Serine";
	codon[59]= "AGA";    aminoacid[59]= "Arginine";
	codon[60]= "AGG";    aminoacid[60]= "Arginine";
	codon[61]= "GGT";    aminoacid[61]= "Glycine";
	codon[62]= "GGC";    aminoacid[62]= "Glycine";
	codon[63]= "GGA";    aminoacid[63]= "Glycine";
	codon[64]= "GGG";    aminoacid[64]= "Glycine";


	//Initialize nucleotides table
	forwardCodonTable[0][0]='A';	forwardCodonTable[0][1]='A';	forwardCodonTable[0][2]='$';	forwardCodonTable[0][3]='$';
	forwardCodonTable[1][0]='C';	forwardCodonTable[1][1]='C';	forwardCodonTable[1][2]='$';	forwardCodonTable[1][3]='$';
	forwardCodonTable[2][0]='G';	forwardCodonTable[2][1]='G';	forwardCodonTable[2][2]='$';	forwardCodonTable[2][3]='$';
	forwardCodonTable[3][0]='T';	forwardCodonTable[3][1]='T';	forwardCodonTable[3][2]='$';	forwardCodonTable[3][3]='$';
	forwardCodonTable[4][0]='R';	forwardCodonTable[4][1]='A';	forwardCodonTable[4][2]='G';	forwardCodonTable[4][3]='$';
	forwardCodonTable[5][0]='Y';	forwardCodonTable[5][1]='C';	forwardCodonTable[5][2]='T';	forwardCodonTable[5][3]='$';
	forwardCodonTable[6][0]='S';	forwardCodonTable[6][1]='G';	forwardCodonTable[6][2]='C';	forwardCodonTable[6][3]='$';
	forwardCodonTable[7][0]='W';	forwardCodonTable[7][1]='A';	forwardCodonTable[7][2]='T';	forwardCodonTable[7][3]='$';
	forwardCodonTable[8][0]='K';	forwardCodonTable[8][1]='G';	forwardCodonTable[8][2]='T';	forwardCodonTable[8][3]='$';
	forwardCodonTable[9][0]='M';	forwardCodonTable[9][1]='A';	forwardCodonTable[9][2]='C';	forwardCodonTable[9][3]='$';
	forwardCodonTable[10][0]='B';	forwardCodonTable[10][1]='C';	forwardCodonTable[10][2]='G';	forwardCodonTable[10][3]='T';
	forwardCodonTable[11][0]='D';	forwardCodonTable[11][1]='A';	forwardCodonTable[11][2]='G';	forwardCodonTable[11][3]='T';
	forwardCodonTable[12][0]='H';	forwardCodonTable[12][1]='A';	forwardCodonTable[12][2]='C';	forwardCodonTable[12][3]='T';
	forwardCodonTable[13][0]='V';	forwardCodonTable[13][1]='A';	forwardCodonTable[13][2]='C';	forwardCodonTable[13][3]='G';
	forwardCodonTable[14][0]='N';	forwardCodonTable[14][1]='$';	forwardCodonTable[14][2]='$';	forwardCodonTable[14][3]='$';

	reverseCodonTable[0][0]='A';	reverseCodonTable[0][1]='T';	reverseCodonTable[0][2]='$';	reverseCodonTable[0][3]='$';
	reverseCodonTable[1][0]='C';	reverseCodonTable[1][1]='G';	reverseCodonTable[1][2]='$';	reverseCodonTable[1][3]='$';
	reverseCodonTable[2][0]='G';	reverseCodonTable[2][1]='C';	reverseCodonTable[2][2]='$';	reverseCodonTable[2][3]='$';
	reverseCodonTable[3][0]='T';	reverseCodonTable[3][1]='A';	reverseCodonTable[3][2]='$';	reverseCodonTable[3][3]='$';
	reverseCodonTable[4][0]='R';	reverseCodonTable[4][1]='T';	reverseCodonTable[4][2]='C';	reverseCodonTable[4][3]='$';
	reverseCodonTable[5][0]='Y';	reverseCodonTable[5][1]='G';	reverseCodonTable[5][2]='A';	reverseCodonTable[5][3]='$';
	reverseCodonTable[6][0]='S';	reverseCodonTable[6][1]='C';	reverseCodonTable[6][2]='G';	reverseCodonTable[6][3]='$';
	reverseCodonTable[7][0]='W';	reverseCodonTable[7][1]='T';	reverseCodonTable[7][2]='A';	reverseCodonTable[7][3]='$';
	reverseCodonTable[8][0]='K';	reverseCodonTable[8][1]='C';	reverseCodonTable[8][2]='A';	reverseCodonTable[8][3]='$';
	reverseCodonTable[9][0]='M';	reverseCodonTable[9][1]='T';	reverseCodonTable[9][2]='G';	reverseCodonTable[9][3]='$';
	reverseCodonTable[10][0]='B';	reverseCodonTable[10][1]='G';	reverseCodonTable[10][2]='C';	reverseCodonTable[10][3]='A';
	reverseCodonTable[11][0]='D';	reverseCodonTable[11][1]='T';	reverseCodonTable[11][2]='C';	reverseCodonTable[11][3]='A';
	reverseCodonTable[12][0]='H';	reverseCodonTable[12][1]='T';	reverseCodonTable[12][2]='G';	reverseCodonTable[12][3]='A';
	reverseCodonTable[13][0]='V';	reverseCodonTable[13][1]='T';	reverseCodonTable[13][2]='G';	reverseCodonTable[13][3]='C';
	reverseCodonTable[14][0]='N';	reverseCodonTable[14][1]='$';	reverseCodonTable[14][2]='$';	reverseCodonTable[14][3]='$';
}
