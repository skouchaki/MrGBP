//
//  feature_extraction.h
//  cluster_matagenom_nt2aa
//
//  Created by Samaneh Kouchaki on 23/05/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef feature_extraction_h
#define feature_extraction_h

#include "oned_lbp.h"

#include <math.h>
#include <map>
#include <vector>


class feature_extraction
{
     std::map <float,int> phases;
public:
    
    void extract(std::string method,std::vector <float> *x,std::string xx,int dimx, int dimy,std::vector <std::vector <float>> *result,int coefn, std::vector <std::string> cp)
    {
        //feature extraction, other methods will be added
       if(method.compare("n_LBP")==0)
        {
            oned_lbp lbp;
            vector <float> r;
            vector <float>r1;
            for (int i=2;i<=coefn;i=i+2)
            {
                lbp.sethist(x,dimx, dimy, i+1,&r);
                copy(r.begin(),r.end(),back_inserter(r1));
                r.clear();
            }
        
            for (int i=0;i<cp.size();i++)
            {
                r1.insert(r1.begin(),(std::stof(cp[i])));
            }
           
            result->push_back(r1);
            r1.clear();
        }
    };
};
#endif /* feature_extraction_h */
