#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>



using namespace std;

//Variable declarations
string readsFolder;
string referenceFile;
string outputFolder;
string outputFile;
string numThreads;
string editDistance;
string additionalBwaFlags;
string additionalPileupFlags;
string file1;
string firstRead[1000];
string secondRead[1000];
string firstReadInPair[1000];
string secondReadInPair[1000];
string singleRead[1000];
string samFileName[1000];

ifstream temp;

char command[1000];


double minimumAlnQual;
double minimumBaseQual;


int readsPaired;
int numberOfFirstRead = 0;
int numberOfSecondRead = 0;
int numberOfFilesInPair = 0;
int numberOfFilesInSingle = 0;
int numberOfSamFiles = 0;
int matePresent;
int a,b;

void checkReadsInFolder(void);
void alignReads(void);
void convertToBam(void);
void sortBam(void);
void mergeBam(void);
void indexReference(void);
void pileupReads(void);
void createFolders(void);
void moveFilesToFolders(void);


int main(int argc, char** argv) {

    readsFolder = argv[1];
    referenceFile = argv[2];
    outputFolder = argv[3];
    outputFile = argv[4];
    editDistance = argv[5];
    numThreads = argv[6];
    additionalBwaFlags = argv[7];
    readsPaired = atoi(argv[8]);
    minimumAlnQual = atoi(argv[9]);
    minimumBaseQual = atoi(argv[10]);
    additionalPileupFlags = argv[11];
    

    createFolders();
    indexReference();
    checkReadsInFolder();
    alignReads();
    convertToBam();
    mergeBam();
    sortBam();
    pileupReads();
    moveFilesToFolders();
    
}


void moveFilesToFolders(void)
{
    strcpy(command,"mv ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.sam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/sam/");
    system(command);
    
    strcpy(command,"rm ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.sai");
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.bam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/bam/");
    system(command);
    
    strcpy(command,"mv *_pileup* ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/pileup/");
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,referenceFile.c_str());
    strcat(command,"* ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/indexedReference");
    system(command);
    
    
}

void createFolders(void)
{
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/sam");
    system(command);
    
    strcpy(command,"mkdir  -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/bam");
    system(command);
    
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/pileup");
    system(command);
    
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/indexedReference");
    system(command);
}

void pileupReads(void)
{
    cout << "Indexing reference for pileup " <<endl;
    cout.flush();
    strcpy(command,"./samtools faidx ");
    strcat(command,referenceFile.c_str());
    system(command);
    cout << command << endl;
    
    cout << "Performing pileup" << endl;
    strcpy(command,"./samtools mpileup -f ");
    strcat(command,referenceFile.c_str());
    strcat(command," ");
    strcat(command,outputFile.c_str());
    strcat(command,"_sorted.bam > ");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileup");
    system(command);
    cout << command << endl;
    
    //Correcting pileup
    cout << "Elaborating pileup " << endl;
    strcpy(command,"./pileupCorr ");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileup ");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileupCorr ");
    cout << command << endl;
    system(command);

}

void indexReference(void)
{
    cout << "Indexing reference file....." << endl;
    cout.flush();
    strcpy(command,"./bwa index ");
    strcat(command,referenceFile.c_str());
    cout << command << endl;
    system(command);
    cout << "Reference successfullly indexed...." << endl;
    cout.flush();

}


void checkReadsInFolder(void)
{
    //Check the reads folder
    strcpy(command,"ls ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.fastq >readsFiles");
    system(command);
    temp.open("readsFiles");
    while(!temp.eof())
    {
        getline(temp,file1,'\n');
        if(!temp.eof())
        {
            cout << file1.substr(file1.size()-8) << endl;
            if (file1.substr(file1.size()-8) == "_1.fastq")
            {
                numberOfFirstRead++;
                firstRead[numberOfFirstRead] = file1.substr(0,file1.size()-8);
            }
            if (file1.substr(file1.size()-8)== "_2.fastq")
            {
                numberOfSecondRead++;
                secondRead[numberOfSecondRead] = file1.substr(0,file1.size()-8);
            }
        }
        
    }
    
    
    
    
    //Split available reads between paired end and single end
    for (a=1;a<=numberOfFirstRead;a++)
    {
        matePresent = 0;
        for(b=1;b<=numberOfSecondRead;b++)
        {
            if (firstRead[a] == secondRead[b])
            {
                matePresent = 1;
                break;
            }
        }
        if (matePresent==1)
        {
            numberOfFilesInPair++;
            firstReadInPair[numberOfFilesInPair]=firstRead[a];
            secondReadInPair[numberOfFilesInPair]=secondRead[b];
        }
        else
        {
            numberOfFilesInSingle++;
            singleRead[numberOfFilesInSingle]= firstRead[a];
        }
        
    }
    
    temp.close();
}

void sortBam(void)
{
    cout << "Sorting bam file " << endl;
    cout.flush();
    
    strcpy(command,"./samtools sort ");
    strcat(command,outputFile.c_str());
    strcat(command,".bam ");
    strcat(command,outputFile.c_str());
    strcat(command,"_sorted");
    system(command);
    cout << command << endl;
}

void convertToBam(void)
{
    cout << "Converting sam file to bam format....." << endl;
    for(a=1;a<=numberOfSamFiles;a++)
    {
        strcpy(command,"./samtools view -bS ");
        strcat(command,samFileName[a].c_str());
        strcat(command,".sam >");
        strcat(command,samFileName[a].c_str());
        strcat(command,".bam");
        system(command);
        cout << command << endl;
    }
}

void mergeBam(void)
{
    cout << "Merging bam files " << endl;
    cout.flush();
    
    strcpy(command,"./samtools merge ");
    strcat(command,outputFile.c_str());
    strcat(command,".bam ");
    for(a=1;a<=numberOfSamFiles;a++)
    {
        strcat(command,samFileName[a].c_str());
        strcat(command,".bam ");
    }
    system(command);
    cout << command << endl;
}



void alignReads(void)
{
    
    // PERFORM THE BWA ALN COMMAND
    cout << "Aligning reads " << endl;
    cout.flush();
    for(a=1;a<=numberOfFilesInPair;a++)
    {
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,firstReadInPair[a].c_str());
        strcat(command,"_1.fastq ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        
        strcat(command," > ");
        strcat(command,firstReadInPair[a].c_str());
        strcat(command,"_1.sai");
        system(command);
        cout << command << endl;
        
        
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,secondReadInPair[a].c_str());
        strcat(command,"_2.fastq ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        strcat(command," > ");
        strcat(command,secondReadInPair[a].c_str());
        strcat(command,"_2.sai");
        system(command);
        cout << command << endl;
    }
    
    for(a=1;a<=numberOfFilesInSingle;a++)
    {
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,singleRead[a].c_str());
        strcat(command,"_1.fastq  ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        
        strcat(command," > ");
        strcat(command,singleRead[a].c_str());
        strcat(command,"_1.sai");
        system(command);
        cout << command << endl;

    }
    
    //PERFORM THE BWA SAMPE/SAMSE COMMAND
    if (readsPaired==1)
    {
        for(a=1;a<=numberOfFilesInPair;a++)
        {
            strcpy(command,"./bwa sampe ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sai ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.fastq ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.fastq > ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = firstReadInPair[a];
            system(command);
            cout << command << endl;
        }
        for(a=1;a<=numberOfFilesInSingle;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,singleRead[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = singleRead[a];
            system(command);
            cout << command << endl;
        }
    }
    else
    {
        for(a=1;a<=numberOfFilesInPair;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = firstReadInPair[a] + "_1";
            system(command);
            cout << command << endl;

            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sai ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.fastq >");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = secondReadInPair[a] + "_2";

            system(command);
            cout << command << endl;
        }
        for(a=1;a<=numberOfFilesInSingle;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,singleRead[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = singleRead[a];

            system(command);
            cout << command << endl;
        }
    }

}




