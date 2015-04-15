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


#include <string>


#ifndef MULTIFASTAANALYZER_H_
#define MULTIFASTAANALYZER_H_

using namespace std;



class multiFastaAnalyzer {


public:
    multiFastaAnalyzer();
    void collectSequences(string);
    void setFileName(string);
    string getSequenceInRange(string, int, int,int);
    string getReverseSequence(string);
    int checkCDSForStartCodon(string);
    int checkCDSForStopCodon(string);
    int checkCDSForIntegerCodons(string);


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
