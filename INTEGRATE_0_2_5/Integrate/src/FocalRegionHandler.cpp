/*
 * FocalRegionHandler.cpp
 *
 *  Created on: Aug 5, 2013
 *      Author: jinzhang
 */

#include "FocalRegionHandler.h"



bool func_sort (region_to_map_t i,region_to_map_t j)
{
	if(i.strand<j.strand)
	{
		return true;
	}
	else if(i.strand==j.strand)
	{
		if(i.tid<j.tid)
		{
			return true;
		}
		else if(i.tid==j.tid)
		{
			if(i.lpos<j.lpos)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}


int FocalRegionHandler::getUion(vector<region_to_map_t> &vtp, vector<region_to_map_t> &vtup) {

	//cout<<"in union"<<endl;
/*
	for(int i=0;i<vtp.size();i++)
	{
		cout<<vtp[i].strand<<" "<<vtp[i].tid<<" "<<vtp[i].lpos<<" "<<vtp[i].rpos<<endl;
	}
*/
	vtup.clear();
	sort(vtp.begin(),vtp.end(),func_sort);
/*
	cout<<"sort"<<endl;
	for(int i=0;i<vtp.size();i++)
	{
		cout<<vtp[i].strand<<" "<<vtp[i].tid<<" "<<vtp[i].lpos<<" "<<vtp[i].rpos<<endl;
	}
*/
	if(vtp.size()==0)
		return 0;
	if(vtp.size()==1)
	{
		vtup.push_back(vtp[0]);
		return 0;
	}
	vtup.push_back(vtp[0]);
	for(int i=1;i<vtp.size();i++)
	{
		if(vtp[i].strand!=vtp[i-1].strand)
		{
			vtup.push_back(vtp[i]);
			continue;
		}
		if(vtp[i].tid!=vtup[vtup.size()-1].tid)
		{
			vtup.push_back(vtp[i]);
			continue;
		}
		if(vtp[i].lpos<=vtup[vtup.size()-1].rpos)
		{
			if(vtp[i].rpos<=vtp[vtup.size()-1].rpos)
			{
				continue;
			}
			else if(vtp[i].rpos>vtp[vtup.size()-1].rpos)
			{
				vtup[vtup.size()-1].rpos=vtp[i].rpos;
				continue;
			}
		}
		else
		{
			vtup.push_back(vtp[i]);
		}
	}
/*
	cout<<"vtup"<<endl;
	for(int i=0;i<vtup.size();i++)
	{
		cout<<vtup[i].strand<<" "<<vtup[i].tid<<" "<<vtup[i].lpos<<" "<<vtup[i].rpos<<endl;
	}
*/
	return 0;
}


