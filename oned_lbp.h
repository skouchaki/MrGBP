//
//  oned_lbp.h
//  cluster_metagenom
//
//  Created by Samaneh Kouchaki on 13/04/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef oned_lbp_h
#define oned_lbp_h

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

class oned_lbp
{
public:
    //calculate 1dLBP of a specific window length
    void sethist(vector <float> *data,int x, int y, int segsize,vector <float>* hist)
    {
        hist->clear();
        for (int i=0;i<(int)pow(2,segsize-1);i++)
            hist->push_back(0);
        float tmp[segsize];

        for (int j=0;j<y;j++)
        {
            for (int i=0;i<x-segsize;i++)
            {
                for (int k=0; k<segsize;k++)
                    tmp[k]=data->at(j*x+k+i);
                int h=cal_hist(tmp,segsize);
                hist->at(h)++;
            }
        }

        for (int i=0;i<pow(2,segsize-1);i++)
            hist->at(i)=(hist->at(i)/((x-segsize)*y));
    }
    int cal_hist(float tmp[],int segsize)
    {
        int s=0;
        int cen=ceil(segsize/2.0);
        for (int i=0;i<cen-1;i++)
        {
            s=s+(1-signbit(tmp[i]-tmp[cen-1]))*pow(2,i)+(1-signbit(tmp[i+cen]-tmp[cen-1]))*pow(2,i+cen-1);
        }
        return s;
    }
};

#endif /* oned_lbp_h */
