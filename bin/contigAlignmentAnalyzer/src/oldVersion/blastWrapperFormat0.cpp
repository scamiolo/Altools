/*
 * blastWrapper.cpp
 *
 *  Created on: Oct 20, 2011
 *      Author: salvo
 */

#include "blastWrapperFormat0.h"
#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include "stdlib.h"

using namespace std;

ifstream inputFile;
blastData myData;
int scannedAlignment;
int partialScannedAlignment;
int endOfFile = 0;
int extraCharacters;

string blastType;

blastWrapper::blastWrapper() {
	// TODO Auto-generated constructor stub

}

void blastWrapper::readHeader(void)
{
	string str;

	getline(inputFile,str,'=');
	inputFile.seekg(-6,ios::cur);
	return;
}

void blastWrapper::convertToLastz(string filename)
{
	int a;
	openFile(filename);
	ofstream outputFile((filename+".cigar").c_str());
	readHeader();
	setBlastType("blastn");


	while(returnEOF()!=1)
	{

		getAlignment();
		cout << "Number of alignment "  <<myData.numberOfAlignments << endl;
		if(returnEOF()!=1)
		{
			for(a=1;a<=myData.numberOfAlignments;a++)
			outputFile << "cigar: " << myData.query << " " << myData.queryStart[a] << " " << myData.queryEnd[a] << " "
					<< myData.queryStrand[a] << " " << myData.alignmentInSubject[a] << " " << myData.subjectStart[a] << " "
					<< myData.subjectEnd[a] << " " << myData.subjectStrand[a] << " " << myData.evalue[a] << " "
					<< myData.cigarCode[a] << endl;
		}


	}
}

string blastWrapper::generateCigarCode(string querySeq, string subjectSeq, string queryStrand, string subjectStrand)
{
	string cigarCode;
	string pos;

	int position;
	int getBackFromInsertion;
	int getBackFromDeletion;
	int totalInsertion;


	//Analyze the case of alignment in the same direction
	if(queryStrand==subjectStrand)
	{
		totalInsertion=0;
			for(position = 1 ; position < subjectSeq.length() ; position++)
			{
				if (char(toupper(querySeq.at(position-1)))!=char(toupper(subjectSeq.at(position-1))) && querySeq.at(position-1)!='-' && subjectSeq.at(position-1) != '-')
				{

					ostringstream oss;
					oss << (position-totalInsertion-1) ;
					pos = oss.str();
					cigarCode += " X " + pos;
				}
				if (querySeq.at(position-1)=='-')
				{
					ostringstream oss;
					oss << (position-totalInsertion-2) ;
					pos = oss.str();
					while(querySeq.at(position-1)=='-')	position++;
					cigarCode += " D " + pos;
				}

				if (subjectSeq.at(position-1) == '-')
				{

					ostringstream oss;
					oss << (position-totalInsertion-2) ;
					pos = oss.str();
					while(subjectSeq.at(position-1)=='-')
					{
						position++;
						totalInsertion++;
					}
					cigarCode += " I " + pos;
				}
			}
			cigarCode.insert(0,"N 0");
	}
	else  //Analyze the case of alignment in different orientation
	{
		totalInsertion=0;
			for(position = subjectSeq.length()-1; position > 0; position--)
			{
				if (char(toupper(querySeq.at(position-1)))!=char(toupper(subjectSeq.at(position-1))) && querySeq.at(position-1)!='-' && subjectSeq.at(position-1) != '-')
				{

					ostringstream oss;
					oss << (subjectSeq.length() - position - totalInsertion) ;
					pos = oss.str();
					cigarCode += " X " + pos;
				}
				if (querySeq.at(position-1)=='-')
				{
					ostringstream oss;
					oss << (subjectSeq.length() -position - totalInsertion-1) ;
					pos = oss.str();
					while(querySeq.at(position-1)=='-')	position--;
					cigarCode += " D " + pos;
				}

				if (subjectSeq.at(position-1) == '-')
				{

					ostringstream oss;
					oss << (subjectSeq.length() - position - totalInsertion-1) ;
					pos = oss.str();
					while(subjectSeq.at(position-1)=='-')
					{
						position--;
						totalInsertion++;
					}
					cigarCode += " I " + pos;
				}
			}
			cigarCode.insert(0,"N 0");
	}



	return cigarCode;
}

void blastWrapper::getAlignmentFormat6(void)
{
	string previousQuery = "null";

	myData.numberOfAlignments = 0;
	while(previousQuery!=myData.query)
	{
		myData.numberOfAlignments++;
		inputFile >> myData.query >> myData.alignmentInSubject[myData.numberOfAlignments] >> myData.identity[myData.numberOfAlignments]
		          >> myData.alignmentLength[myData.numberOfAlignments] >> myData.numberOfMismatches[myData.numberOfAlignments]
		          >> myData.numberOfGaps[100000] >> myData.queryStart[myData.numberOfAlignments] >> myData.queryEnd[myData.numberOfAlignments]
		          >> myData.subjectStart[myData.numberOfAlignments] >> myData.subjectEnd[myData.numberOfAlignments]
		          >> myData.evalue[myData.numberOfAlignments] >> myData.score[myData.numberOfAlignments];

		previousQuery = myData.query;


	}

}

void blastWrapper::getAlignment(void)
{
	int positionInFile;

	string str;
	char c;
	int a;
	int subjectStartPositionObtained;




	//start scanning

	//get query name
	getline(inputFile,str,' ');

	if(str==":" || str=="" ) //end of File reached
	{
		endOfFile=1;
		return;
	}
	getline(inputFile, myData.query, ' ');
	getline(inputFile,str,'\n');
	if (inputFile.eof()) return;

	cout << "query name : " << myData.query << endl;

	//get query length
	c='1';
	while(c!='=') c = inputFile.get();
	getline(inputFile,str,'\n');
	myData.queryLength = atoi(str.c_str());

	//get alignment information
	myData.numberOfAlignments=0;

	//alignment name
	while(c!='>') c = inputFile.get();
	c = inputFile.get();
	begin:
	
	myData.numberOfAlignments++;
	getline(inputFile, myData.alignmentInSubject[myData.numberOfAlignments],'\n');
	//cout << myData.alignmentInSubject[myData.numberOfAlignments] << endl;

	//suject length
	moveToDoubleNewLine();
	getline(inputFile,str);
	myData.subjectLength[myData.numberOfAlignments] = atoi(str.c_str());
	//cout << "Sequence length " << myData.subjectLength[myData.numberOfAlignments] << endl;
	//alignment statistics
	//Score
	score:
	c='1';
	while(c!='=') c = inputFile.get();
	c=' ';
	while(c==' ') c = inputFile.get();
	inputFile.seekg(-1,ios::cur);

	getline(inputFile,str,' ');
	myData.score[myData.numberOfAlignments] = atof(str.c_str());
	//cout << "found score " << myData.score[myData.numberOfAlignments] << endl;


	//bits
	getline(inputFile,str,'(');
	getline(inputFile,str,')');
	myData.bits[myData.numberOfAlignments] = atoi(str.c_str());
	//evalue
	getline(inputFile,str,'=');
	c = inputFile.get();
	getline(inputFile,str,'\n');
	myData.evalue[myData.numberOfAlignments] = atof(str.c_str());
	//identity
	getline(inputFile,str,'=');
	c = inputFile.get();
	getline(inputFile,str,'/');
	myData.identity[myData.numberOfAlignments] = atoi(str.c_str());
	getline(inputFile,str,' ');
	myData.alignmentLength[myData.numberOfAlignments] = atoi(str.c_str());
	//positive (if any)
	if(blastType=="blastx")
	{
		getline(inputFile,str,'=');
		c = inputFile.get();
		getline(inputFile,str,'/');
		myData.positive[myData.numberOfAlignments] = atoi(str.c_str());
	}

	//gaps
	getline(inputFile,str,'=');
	c = inputFile.get();
	getline(inputFile,str,'/');
	myData.gaps[myData.numberOfAlignments] = atoi(str.c_str());
	getline(inputFile,str,'\n');

	//frame (if any)
	if(blastType=="blastx")
	{
		getline(inputFile,str,'=');
		myData.frame[myData.numberOfAlignments] = atoi(str.c_str());
		while(c!='y') c = inputFile.get();
	}

	//strnad (if any)
	if(blastType=="blastn")
	{
		getline(inputFile,str,'=');
		getline(inputFile,myData.strand[myData.numberOfAlignments],'\n');
		while(c!='y') c = inputFile.get();
	}


	//alignment
	subjectStartPositionObtained=0;
	c = inputFile.get();
	c = inputFile.get();
	getline(inputFile,str,' ');
	myData.queryStart[myData.numberOfAlignments] = atoi(str.c_str());

	//cout << "query start : " << myData.queryStart[myData.numberOfAlignments] << endl;
	//cout << "alignment length : " << myData.alignmentLength[myData.numberOfAlignments] << endl;
	c=' ';
	while(c==' ') c = inputFile.get();
	inputFile.seekg(-1,ios::cur);

	scannedAlignment = 0;
	partialScannedAlignment=0;
	myData.querySequence[myData.numberOfAlignments] = "";
	myData.subjectSequence[myData.numberOfAlignments] = "";
	myData.consensusSequence[myData.numberOfAlignments] = "";

	while(scannedAlignment< myData.alignmentLength[myData.numberOfAlignments])
	{
		//get the query sequence

		c = inputFile.get();

		/*if (partialScannedAlignment==0)
		{
			cout << c ;
			getchar();
		}*/



		if (c!= ' ')
		{
			myData.querySequence[myData.numberOfAlignments] += c;
			scannedAlignment++;
			partialScannedAlignment++;
		}

		if(c==' ') //get the consensus sequence
		{
			extraCharacters=0;
			getline(inputFile,str,'\n');
			c=' ';
			while(c==' ') c = inputFile.get();
			myData.consensusSequence[myData.numberOfAlignments] += c;
			for(a=0;a<partialScannedAlignment-1;a++)
			{
				c = inputFile.get();
				if(c!='\n') myData.consensusSequence[myData.numberOfAlignments] += c;
				else
				{
					myData.consensusSequence[myData.numberOfAlignments].insert(myData.consensusSequence[myData.numberOfAlignments].length()-partialScannedAlignment+1," ");
					extraCharacters++;
				}
				inputFile.seekg(-extraCharacters,ios::cur);
			}
			//get the subject sequence
			getline(inputFile,str, ' ');
			c = inputFile.get();
			getline(inputFile,str, ' ');
			if(subjectStartPositionObtained == 0)
			{
				myData.subjectStart[myData.numberOfAlignments] = atoi(str.c_str());
				subjectStartPositionObtained = 1;
			}
			c=' ';
			while (c==' ') c = inputFile.get();
			myData.subjectSequence[myData.numberOfAlignments] += c;
			for(a=1;a<partialScannedAlignment;a++) myData.subjectSequence[myData.numberOfAlignments] += inputFile.get();
			partialScannedAlignment = 0;
			getline(inputFile,str, 'y');
			c=' ';
			while (c==' ') c = inputFile.get();
			getline(inputFile,str, ' ');
			c=' ';
			while (c==' ') c = inputFile.get();
			inputFile.seekg(-1,ios::cur);
			//getchar();
		}

		if(scannedAlignment == myData.alignmentLength[myData.numberOfAlignments])
		{
			//get last bit of consensus sequence
			getline(inputFile,str, ' ');
			getline(inputFile,str, ' ');
			getline(inputFile,str, ' ');
			extraCharacters=0;
			myData.queryEnd[myData.numberOfAlignments] = atoi(str.c_str());
			c=' ';
			while(c==' ') c = inputFile.get();
			myData.consensusSequence[myData.numberOfAlignments] += c;

			for(a=0;a<partialScannedAlignment-1;a++)
			{
				c = inputFile.get();
				if(c!='\n') myData.consensusSequence[myData.numberOfAlignments] += c;
				else
				{
					myData.consensusSequence[myData.numberOfAlignments].insert(myData.consensusSequence[myData.numberOfAlignments].length()-partialScannedAlignment+1," ");
				    extraCharacters++;
				}
				inputFile.seekg(-extraCharacters,ios::cur);
			}

			//get last bit of subject sequence
			getline(inputFile,str, ' ');
			getline(inputFile,str, ' ');
			getline(inputFile,str, ' ');
			c=' ';
			while(c==' ') c = inputFile.get();
			myData.subjectSequence[myData.numberOfAlignments] += c;

			for(a=0;a<partialScannedAlignment;a++) myData.subjectSequence[myData.numberOfAlignments] += inputFile.get();
			c = inputFile.get();
			getline(inputFile,str, '\n');

			myData.subjectEnd[myData.numberOfAlignments] = atoi(str.c_str());
		}


	}
	//cout << endl;


	//cout << myData.querySequence[myData.numberOfAlignments] << endl;
	//cout << myData.consensusSequence[myData.numberOfAlignments] << endl;
	//cout << myData.subjectSequence[myData.numberOfAlignments] << endl;
	//cout << "Query: " << myData.queryStart[myData.numberOfAlignments] << "-" << myData.queryEnd[myData.numberOfAlignments] << endl;
	//cout << "Subject: " << myData.subjectStart[myData.numberOfAlignments] << "-" << myData.subjectEnd[myData.numberOfAlignments] << endl;


	//define strands
	if(myData.queryStart[myData.numberOfAlignments] > myData.queryEnd[myData.numberOfAlignments])
	     myData.queryStrand[myData.numberOfAlignments] = "-";
	else myData.queryStrand[myData.numberOfAlignments] = "+";

	if(myData.subjectStart[myData.numberOfAlignments] > myData.subjectEnd[myData.numberOfAlignments])
	     myData.subjectStrand[myData.numberOfAlignments] = "-";
	else myData.subjectStrand[myData.numberOfAlignments] = "+";

	myData.cigarCode[myData.numberOfAlignments] = generateCigarCode(myData.querySequence[myData.numberOfAlignments],
	                                                                myData.subjectSequence[myData.numberOfAlignments],
	                                                                myData.queryStrand[myData.numberOfAlignments],
	                                                                myData.subjectStrand[myData.numberOfAlignments]);
	for(a=0;a<3;a++)
	{
		c = inputFile.get();
		if (inputFile.eof())
		{
			endOfFile = 1;
			break;
		}
	}
	if (c=='>')
	{
		c = inputFile.get();
		if (inputFile.eof())
		{
			endOfFile = 1;
			return;
		}
		goto begin;
	}
	else
	{
		c = inputFile.get();
		if (c=='S')
		{
			myData.numberOfAlignments++;
			myData.alignmentInSubject[myData.numberOfAlignments] = myData.alignmentInSubject[myData.numberOfAlignments-1];
			//getline(inputFile,str,'=');
			//c=inputFile.get();
			//c=inputFile.get();
			goto score;
		}
		else readHeader();
	}

}

void blastWrapper::setBlastType(string type)
{
	blastType = type;
	return;
}

int blastWrapper::returnEOF(void)
{
	return endOfFile;
}
void blastWrapper::moveToDoubleNewLine(void)
{
	char c;
	int flag = 0;


	c = inputFile.get();
	while(c!='\n' && flag==0)
	{
		c = inputFile.get();
		if(c=='\n')
		{
			c = inputFile.get();
			if(c=='\n') flag = 1;
		}
	}

	while(c!='=')
	{
		inputFile.seekg(-2,ios::cur);
		c = inputFile.get();
	}

}

void blastWrapper::openFile(string filename)
{
	inputFile.open(filename.c_str());
	return;
}

blastWrapper::~blastWrapper() {
	// TODO Auto-generated destructor stub
}

