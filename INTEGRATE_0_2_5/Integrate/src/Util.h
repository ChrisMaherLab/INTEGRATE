/*
 * Util.h
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <vector>
#include <iostream>

using namespace std;


uint32_t getFilelength(char *file);
int readBlock(char * block, int length, FILE *infile);


extern map<int,char> intChar;
extern map<char,char> charChar;

extern map<string,char> tableAmino;


int InitialIntChar();


char getCharComp(char reada);
char getCharA(int reada);

int getPeptide(vector<char> & seq5p, vector<char> & seq, int start_pos, vector<char> & peptide, int & full, int & left);

#endif /* UTIL_H_ */
