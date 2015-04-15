/* Copyright 2013 Salvatore Camiolo
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*/

#include "multiFastaAnalyzer.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>




multiFastaAnalyzer::multiFastaAnalyzer() {
    // TODO Auto-generated constructor stub

    }

void multiFastaAnalyzer::collectSequences(string fileName) //load all the sequences in memory together with thei position in the file and the locus name
{
    string preSequence;
    int a;
    sequence = new string[500000];
    locus = new string[500000];
    title = new string[500000];
    positionInFile = new int[500000];


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

            if (fileToOpen.eof())
            {
                fileToOpen.seekg(-1,ios::cur);
                for(a=0;a<preSequence.length();a++)
                {
                    if (preSequence[a]!='\n') sequence[numberOfSequences] += toupper(preSequence[a]);
                }
                break;
            }

            fileToOpen.seekg(-1,ios::cur);
            for(a=0;a<preSequence.length();a++)
            {
                if (preSequence[a]!='\n') sequence[numberOfSequences] += toupper(preSequence[a]);
            }

            if (fmod(numberOfSequences,1000)==0) cout << numberOfSequences << " loaded from file " << fileName << endl;
        }
        cout << numberOfSequences << " loaded from file " << fileName << endl;
    }
    else {
            cout << "The file " << fileName << "does not exist!\nThe program will now exit.\n" << endl;
            getchar();
            exit(0);
        }


    fileToOpen.close();

return;
}


string multiFastaAnalyzer::getSequenceInRange(string chromosome, int start, int end,int strand)
{
    string chromosomeSequence;
    string searchedSequence;
    string reverseSearchedSequence;


    chromosomeSequence = retrieveSequence(chromosome);


    if( (start<chromosomeSequence.length()) && (end <= chromosomeSequence.length()) && chromosome !="")
    {
        if(strand==0)
        {
            searchedSequence = chromosomeSequence.substr(start-1,end-start+1);
            return searchedSequence;
        }
        else
        {
            searchedSequence = chromosomeSequence.substr(start-1,end-start+1);
            reverseSearchedSequence = getReverseSequence(searchedSequence);
            return reverseSearchedSequence;
        }
    }
    else return "Sequence_is_Out_Of_Chromosome";
}

string multiFastaAnalyzer::getReverseSequence(string sequenceToReverse)
{
    string reversedSequence;
    int a,i;

    i=0;
    for(a=sequenceToReverse.length()-1;a>=0;a--)
    {
        if(sequenceToReverse.at(a)=='A' || sequenceToReverse.at(a)=='G' || sequenceToReverse.at(a)=='C' || sequenceToReverse.at(a)=='T' ||
                sequenceToReverse.at(a)=='a' || sequenceToReverse.at(a)=='g' || sequenceToReverse.at(a)=='c' || sequenceToReverse.at(a)=='t')
        {
            if(sequenceToReverse.at(a)=='A' || sequenceToReverse.at(a)=='a') reversedSequence.append(1,'T');
            if(sequenceToReverse.at(a)=='G' || sequenceToReverse.at(a)=='g') reversedSequence.append(1,'C');
            if(sequenceToReverse.at(a)=='T' || sequenceToReverse.at(a)=='t') reversedSequence.append(1,'A');
            if(sequenceToReverse.at(a)=='C' || sequenceToReverse.at(a)=='c') reversedSequence.append(1,'G');
            i++;
        }
        else
        {
            reversedSequence.append(1,sequenceToReverse.at(a));
            i++;
        }
    }

    return reversedSequence;
}


int multiFastaAnalyzer::checkCDSForStartCodon(string cds)
{
    string codon;
    if (cds.length()<3) return 1;
    codon.append(1,cds.at(0));
    codon.append(1,cds.at(1));
    codon.append(1,cds.at(2));

    if (codon=="ATG") return 0;
    else return 1;
}

int multiFastaAnalyzer::checkCDSForStopCodon(string cds)
{
    string codon;
    if (cds.length()<3) return 1;
    codon.append(1,cds.at(cds.length()-3));
    codon.append(1,cds.at(cds.length()-2));
    codon.append(1,cds.at(cds.length()-1));

    if (codon=="TGA" || codon=="TAG" || codon=="TAA") return 0;
    else return 1;
}

int multiFastaAnalyzer::checkCDSForIntegerCodons(string cds)
{
    int reminder;
    reminder = fmod(cds.length(),3);

    if (reminder == 0) return 0;
    else return 1;
}


string multiFastaAnalyzer::retrieveSequence(string locusToRetrieve)
{
    int i;

    for (i=0;i<=numberOfSequences;i++)
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
