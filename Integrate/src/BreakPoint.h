/*
 * BreakPoint.h
 *
 *  Created on: March 11, 2014
 *      Author: jinzhang
 */

#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_

#include <iostream>
#include <stdint.h>
#include <fstream>
#include <ctime>
#include <sstream>
#include <algorithm>

using namespace std;


#include "MyTypes.h"
#include "Alignment.h"
#include "Reference.h"



class BreakPoint {
public:

    
    BreakPoint();
    virtual ~BreakPoint();
	
    int getBreakPoints(vector<break_point_record_t> & bkvec, char * filename,char * filename1,char * filename2,char * refname, Reference & ref, char* sample_name);
    
    int getOneBKRNA(break_point_record_t & bkt,Reference & ref);
  
    int getOneBKDNA(break_point_record_t & bkt,Reference & ref);

 
    int printOneBK(break_point_record_t & bkt, Reference & ref, ofstream & outFile);
    int AddOneVCF(break_point_record_t & bkt,  Reference & ref);
    int printOneBEDPE(break_point_record_t & bkt, Reference & ref, ofstream & outFile);
    
};

#endif /* BREAKPOINT_H_ */
