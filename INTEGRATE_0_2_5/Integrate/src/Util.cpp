/*
 * Util.cpp
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 */

#include "Util.h"


/*
 * count the length of chars in a file
 *
 *
 */


uint32_t getFilelength(char *file)
{
	struct stat filestatus;
    stat( file, &filestatus );
    uint32_t length = filestatus.st_size;
    return length;
}


/*
 * Read a block of chars from file till a \n of file
 * return the actual length;
 *
 */

int readBlock(char * block, int length, FILE *infile)
{
	fread(block,1,length,infile);
	if(block[length-1]!='\n' && block[length-1]!=EOF)
	{
		fgets(block+length,1024,infile);
		length+=strlen(block+length);
	}


	return length;
}


int InitialIntChar()
{
	intChar.insert(pair<int,char>(1,'A'));
	intChar.insert(pair<int,char>(2,'C'));
	intChar.insert(pair<int,char>(4,'G'));
	intChar.insert(pair<int,char>(8,'T'));
	intChar.insert(pair<int,char>(15,'N'));

	charChar.insert(pair<char,char>('A','T'));
	charChar.insert(pair<char,char>('C','G'));
	charChar.insert(pair<char,char>('G','C'));
	charChar.insert(pair<char,char>('T','A'));
	charChar.insert(pair<char,char>('N','N'));
    
    
    tableAmino.insert(pair<string,char>("GCT",'A'));
    tableAmino.insert(pair<string,char>("GCC",'A'));
    tableAmino.insert(pair<string,char>("GCA",'A'));
    tableAmino.insert(pair<string,char>("GCG",'A'));
    tableAmino.insert(pair<string,char>("CGT",'R'));
    tableAmino.insert(pair<string,char>("CGC",'R'));
    tableAmino.insert(pair<string,char>("CGA",'R'));
    tableAmino.insert(pair<string,char>("CGG",'R'));
    tableAmino.insert(pair<string,char>("AGA",'R'));
    tableAmino.insert(pair<string,char>("AGG",'R'));
    tableAmino.insert(pair<string,char>("AAT",'N'));
    tableAmino.insert(pair<string,char>("AAC",'N'));
    tableAmino.insert(pair<string,char>("GAT",'D'));
    tableAmino.insert(pair<string,char>("GAC",'D'));
    tableAmino.insert(pair<string,char>("TGT",'C'));
    tableAmino.insert(pair<string,char>("TGC",'C'));
    tableAmino.insert(pair<string,char>("CAA",'Q'));
    tableAmino.insert(pair<string,char>("CAG",'Q'));
    tableAmino.insert(pair<string,char>("GAA",'E'));
    tableAmino.insert(pair<string,char>("GAG",'E'));
    tableAmino.insert(pair<string,char>("GGT",'G'));
    tableAmino.insert(pair<string,char>("GGC",'G'));
    tableAmino.insert(pair<string,char>("GGA",'G'));
    tableAmino.insert(pair<string,char>("GGG",'G'));
    tableAmino.insert(pair<string,char>("CAT",'H'));
    tableAmino.insert(pair<string,char>("CAC",'H'));
    tableAmino.insert(pair<string,char>("ATT",'I'));
    tableAmino.insert(pair<string,char>("ATC",'I'));
    tableAmino.insert(pair<string,char>("ATA",'I'));
    tableAmino.insert(pair<string,char>("TTA",'L'));
    tableAmino.insert(pair<string,char>("TTG",'L'));
    tableAmino.insert(pair<string,char>("CTT",'L'));
    tableAmino.insert(pair<string,char>("CTC",'L'));
    tableAmino.insert(pair<string,char>("CTA",'L'));
    tableAmino.insert(pair<string,char>("CTG",'L'));
    tableAmino.insert(pair<string,char>("AAA",'K'));
    tableAmino.insert(pair<string,char>("AAG",'K'));
    tableAmino.insert(pair<string,char>("ATG",'M'));
    tableAmino.insert(pair<string,char>("TTT",'F'));
    tableAmino.insert(pair<string,char>("TTC",'F'));
    tableAmino.insert(pair<string,char>("CCT",'P'));
    tableAmino.insert(pair<string,char>("CCC",'P'));
    tableAmino.insert(pair<string,char>("CCA",'P'));
    tableAmino.insert(pair<string,char>("CCG",'P'));
    tableAmino.insert(pair<string,char>("TCT",'S'));
    tableAmino.insert(pair<string,char>("TCC",'S'));
    tableAmino.insert(pair<string,char>("TCA",'S'));
    tableAmino.insert(pair<string,char>("TCG",'S'));
    tableAmino.insert(pair<string,char>("AGT",'S'));
    tableAmino.insert(pair<string,char>("AGC",'S'));
    tableAmino.insert(pair<string,char>("ACT",'T'));
    tableAmino.insert(pair<string,char>("ACC",'T'));
    tableAmino.insert(pair<string,char>("ACA",'T'));
    tableAmino.insert(pair<string,char>("ACG",'T'));
    tableAmino.insert(pair<string,char>("TGG",'W'));
    tableAmino.insert(pair<string,char>("TAT",'Y'));
    tableAmino.insert(pair<string,char>("TAC",'Y'));
    tableAmino.insert(pair<string,char>("GTT",'V'));
    tableAmino.insert(pair<string,char>("GTC",'V'));
    tableAmino.insert(pair<string,char>("GTA",'V'));
    tableAmino.insert(pair<string,char>("GTG",'V'));
    tableAmino.insert(pair<string,char>("TAA",'X'));//X for stop
    tableAmino.insert(pair<string,char>("TGA",'X'));
    tableAmino.insert(pair<string,char>("TAG",'X'));
    
	return 0;
}

char getCharComp(char reada)
{
	return charChar[reada];

}



char getCharA(int reada)
{
	return intChar[reada];
}












char getAmino(string a)
{
    return tableAmino[a];
}






int getPeptide(vector<char> & seq5p, vector<char> & seq, int start_pos, vector<char> & peptide, int & full, int & left)
{
//cout<<"in getPeptide "<<seq5p.size()<<" "<<seq.size()<<endl;
    full=0;
//cout<<"start_pos "<<start_pos<<endl;
    for(int i=start_pos-1;i<seq5p.size()-2;i+=3)
    {
        full++;
    }
    
    left=seq5p.size()-3*full-start_pos+1;

//cout<<"full "<<full<<" left"<<left<<endl;
 
    for(int i=start_pos-1;i<seq.size()-2;i+=3)
    {
        string a;
        for(int j=0;j<3;j++)
            a=a+seq[i+j];
        char am=getAmino(a);
//cout<<"char "<<am<<endl;
        if(am!='X')
            peptide.push_back(am);
        else
            break;
    }
//cout<<"out"<<endl;
    return 0;
}






