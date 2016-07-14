/*
 * FocalRegionHandler.h
 *
 *  Created on: Aug 2, 2013
 *      Author: jinzhang
 */

#ifndef FOCALREGIONHANDLER_H_
#define FOCALREGIONHANDLER_H_



#include <iostream>
#include "MyTypes.h"
#include <vector>
#include <algorithm>
#include <iterator>


using namespace std;

//Could use interval tree;

class FocalRegionHandler {

public:
	int getUion(vector<region_to_map_t> &vtp, vector<region_to_map_t> &vtup);
};


#endif /* FOCALREGIONHANDLER_H_ */
