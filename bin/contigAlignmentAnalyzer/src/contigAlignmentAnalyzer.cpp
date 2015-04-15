//============================================================================
// Name        : contigAlignmentAnalyzer.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "completePileupAnalyzer.h"



using namespace std;

string inputFileName;
string filename;
string outputFileName;
string outputFolderName;
string pileupFolderName1;
string pileupFolderName2;
string alignmentFormat;
string freeSlot1;
string speciesName1;
string speciesName2;

int minimumBitCutOff;
int check_endOfFile;
int minimumInsertionSize;
int minimumDeletionSize;
int minimumTranslocationSize;
int minimumInversionSize;
int minimumTranslocationPlusInversionSize;
int minimumDuplicationSize;
int searchGoldenPolymorphysmSet;
int generateInsertionsGE;
int generateDeletionsGE;
int generateTranslocationGE;
int generateInversionGE;
int generateTranslocPlusInversion;
int generateDuplicationGE;
int foundCommon;
int a,b,c;
int (*numSnp1)[10000000] = new int[100][10000000];
int (*numSnp2)[10000000] = new int[100][10000000];


ofstream uniques1;
ofstream uniques2;
ofstream common;


completePileupAnalyzer myPileup1;
completePileupAnalyzer myPileup2;



int main(int argc, char** argv) {


    pileupFolderName1 = argv[1];
    pileupFolderName2 = argv[2];
    outputFolderName = argv[3];
    speciesName1 = argv[4];
    speciesName2 = argv[5];

    filename = outputFolderName + "/" + speciesName1 + "_uniqPolymorphysms.txt" ;
    uniques1.open(filename.c_str());
    filename = outputFolderName + "/" + speciesName2 + "_uniqPolymorphysms.txt" ;
    uniques2.open(filename.c_str());
    filename = outputFolderName + "/commonPolymorphisms.txt" + speciesName2 + "_uniqPolymorphysms.txt" ;
    common.open(filename.c_str());

    
    myPileup1.setFolder(pileupFolderName1);
    myPileup1.collectPolymorphysms();
    
    myPileup2.setFolder(pileupFolderName2);
    myPileup2.collectPolymorphysms();

    
    for(a=1;a<=myPileup1.numberOfFiles;a++)
    {
        cout << "Comparing species 1 SNPs species 2 SNPs on chromosome " << myPileup1.fileInPileupFolder[a]<< endl;
        //Compare snp species 1 to snps species 2
        for(b=1; b<=myPileup1.numberOfSnp[a];b++)
        {
            foundCommon = 0;
            for(c=1;c<=myPileup2.numberOfSnp[a];c++)
            {
                if(myPileup1.snpPosition[a][b] == myPileup2.snpPosition[a][c])
                {
                    
                    cout << "Found concordant SNP at " <<myPileup1.snpPosition[a][b] << "             \r";
                    common << myPileup1.fileInPileupFolder[a] << "\tSNP\t" << myPileup1.snpPosition[a][b] << endl;
                    foundCommon = 1;
                    break;
                }
            }
        
         if (foundCommon==0) uniques1 << myPileup1.fileInPileupFolder[a] << "\tSNP\t" << myPileup1.snpPosition[a][b] << endl;
        }
        
       
            
            //Compare snp species 2 to snps species 1
            cout << "Comparing species 2 SNPs species 1 SNPs on chromosome "<< myPileup1.fileInPileupFolder[a] << endl;
            for(b=1; b<=myPileup2.numberOfSnp[a];b++)
            {
                foundCommon = 0;

                for(c=1;c<=myPileup1.numberOfSnp[a];c++)
                {
                    if(myPileup2.snpPosition[a][b] == myPileup1.snpPosition[a][c])
                    {
                        
                        cout << "Found concordant snp at " <<myPileup1.snpPosition[a][b] << "         \r";
                        common << myPileup1.fileInPileupFolder[a] << "\tSNP\t" << myPileup1.snpPosition[a][b] << endl;
                        foundCommon = 1;
                        break;
                    }
                }
           
            if (foundCommon==0) uniques2 << myPileup2.fileInPileupFolder[a] << "\tSNP\t" << myPileup2.snpPosition[a][b] << endl;
            }
        
        
        cout << "Comparing species 1 INDELs species 2 INDELs on chromosome " << myPileup1.fileInPileupFolder[a]<< endl;
        //Compare snp species 1 to snps species 2
        for(b=1; b<=myPileup1.numberOfIndel[a];b++)
        {
            foundCommon = 0;
            for(c=1;c<=myPileup2.numberOfIndel[a];c++)
            {
                if(myPileup1.indelPosition[a][b] == myPileup2.indelPosition[a][c])
                {
                    
                    cout << "Found concordant INDEL at " <<myPileup1.indelPosition[a][b] << "             \r";
                    common << myPileup1.fileInPileupFolder[a] << "\tINDEL\t" << myPileup1.indelPosition[a][b] << endl;
                    foundCommon = 1;
                    break;
                }
            }
            
            if (foundCommon==0) uniques1 << myPileup1.fileInPileupFolder[a] << "\tINDEL\t" << myPileup1.indelPosition[a][b] << endl;
        }
        
        
        
        //Compare snp species 2 to snps species 1
        cout << "Comparing species 2 INDELs species 1 INDELs on chromosome "<< myPileup1.fileInPileupFolder[a] << endl;
        for(b=1; b<=myPileup2.numberOfIndel[a];b++)
        {
            foundCommon = 0;
            
            for(c=1;c<=myPileup1.numberOfIndel[a];c++)
            {
                if(myPileup2.indelPosition[a][b] == myPileup1.indelPosition[a][c])
                {
                    
                    cout << "Found concordant INDEL at " <<myPileup1.indelPosition[a][b] << "         \r";
                    common << myPileup1.fileInPileupFolder[a] << "\tINDEL\t" << myPileup1.indelPosition[a][b] << endl;
                    foundCommon = 1;
                    break;
                }
            }
            
            if (foundCommon==0) uniques2 << myPileup2.fileInPileupFolder[a] << "\tINDEL\t" << myPileup2.indelPosition[a][b] << endl;
        }
    }
    
    cout << "The process successfully terminated!!" << endl;

	return 0;
}
