/*
 * Updator.h
 *
 *  Created on: Jun 26, 2016
 *      Author: jinzhang
 */

#ifndef UPDATOR_H_
#define UPDATOR_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <fstream>


using namespace std;


#include "MyTypes.h"
#include "Reference.h"
#include "Gene.h"

class Updator {

public:
    
    Updator(){};
    ~Updator(){};
    
    int update(int gid5p, int gid3p, uint32_t & pos5p, uint32_t & pos3p, split_rna_t & st, Reference & ref, Gene & g);
};    




#endif 
