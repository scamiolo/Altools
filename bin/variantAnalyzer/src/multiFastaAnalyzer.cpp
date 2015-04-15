/*
 * multiFastaAnalyzer.cpp
 *
 *  Created on: Sep 8, 2011
 *      Author: salvo
 */

#include "multiFastaAnalyzer.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>




multiFastaAnalyzer::multiFastaAnalyzer() {
	// TODO Auto-generated constructor stub

	}

void multiFastaAnalyzer::collectSequences(string fileName) //load all the sequences in memory together with thei position in the file and the locus name
{
	string preSequence;
	int a;
	sequence = new string[1000000];
	locus = new string[1000000];
	title = new string[1000000];
	positionInFile = new int[1000000];


	numberOfSequences = -1;
	ifstream fileToOpen(fileName.c_str());

	if(fileToOpen)
	{
		while(!fileToOpen.eof())
		{
			numberOfSequences++;

			//get locus and title of the sequence
			getline(fileToOpen,title[numberOfSequences]);
			locus[numberOfSequences] = title[numberOfSequences].substr(1,title[numberOfSequences].find_first_of(" \t\n")-1);

			//get the sequence string, position in file and remove new lines from the sequence string
			positionInFile[numberOfSequences] = fileToOpen.tellg();
			getline(fileToOpen,preSequence,'>');
			if (fileToOpen.eof()) break;
			fileToOpen.seekg(-1,ios::cur);
			for(a=0;a<preSequence.length();a++)
			{
				if (preSequence[a]!='\n') sequence[numberOfSequences] += preSequence[a];
			}

			if (fmod(numberOfSequences,1000)==0) cout << numberOfSequences << " loaded from file " << fileName << endl;
		}
	}
	else
	{
		cout << "The file " << fileName << "does not exist!\nThe program will now exit.\n" << endl;
		exit(0);
	}

	fileToOpen.close();

return;
}


string multiFastaAnalyzer::retrieveSequence(string locusToRetrieve)
{
	int i;

	for (i=0;i<numberOfSequences;i++)
	{
		if(strcmp(locusToRetrieve.c_str(),locus[i].c_str())==0)
		{
			return sequence[i];
		}
	}
	return "None";
}



void multiFastaAnalyzer::setFileName(string fileName) // set the file name to work with
{
	fileToAnalyze = fileName;
	return;
}

multiFastaAnalyzer::~multiFastaAnalyzer() {
	// TODO Auto-generated destructor stub
}
