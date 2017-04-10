//
//  nucreps.h
//  cluster_matagenom_nt2aa
// numerically represent data
//  Created by Samaneh Kouchaki on 25/05/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef nucreps_h
#define nucreps_h

#include <string>
#include <vector>
#include <map>
#include <math.h>

class nucreps
{
    std::vector <std::map<char,std::vector <float>>> mreps;
    std::map <std::string,float> repsl;
    
public:
    nucreps()
    {
        repsl["Atomic"] = 0;
        repsl["Complex"] = 1;
        repsl["Dna_Walk"] = 2;
        repsl["EIIP"] = 3;
        repsl["Integer"] = 4;
        repsl["Paired"] = 5;
        repsl["Real"] = 6;
        repsl["Tetrahedron"] = 7;
        repsl["Voss"] = 8;
        repsl["Z_curve"] = 9;
        
        std::map<char,std::vector <float>> mrep;
        std::vector <float> t;//=new double[4];
        t.clear();  
        //(0)
        t.push_back(70);
        mrep['A']=t;
        t[0]=58;
        mrep['C']=t;
        t[0]=66;
        mrep['G']=t;
        t[0]=78;
        mrep['T']=t;
        t[0]=0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(1)
        t.clear();
        t.push_back(1.0);
        t.push_back(1.0);
        mrep['A']=t;
        t[0]=-1.0;
        t[1]=1.0;
        mrep['C']=t;
        t[0]=-1.0;
        t[1]=-1.0;
        mrep['G']=t;
        t[0]=1.0;
        t[1]=-1.0;
        mrep['T']=t;
        t[0]=0.0;
        t[1]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(2)
        t.push_back(0);
        t.push_back(-1.0);
        mrep['A']=t;
        t[0]=1.0;
        t[1]=0;
        mrep['C']=t;
        t[0]=-1.0;
        t[1]=0;
        mrep['G']=t;
        t[0]=0;
        t[1]=1.0;
        mrep['T']=t;
        t[0]=0.0;
        t[1]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(3)
        t.clear();
        t.push_back(0.1260);
        mrep['A']=t;
        t[0]=0.1340;
        mrep['C']=t;
        t[0]=0.0806;
        mrep['G']=t;
        t[0]=0.1335;
        mrep['T']=t;
        t[0]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(4)
        t.clear();
        t.push_back(2.0);
        mrep['A']=t;
        t[0]=-1;
        mrep['C']=t;
        t[0]=1;
        mrep['G']=t;
        t[0]=-2;
        mrep['T']=t;
        t[0]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(5)
        t.clear();
        t.push_back(1.0);
        mrep['A']=t;
        t[0]=-1.0;
        mrep['C']=t;
        t[0]=-1.0;
        mrep['G']=t;
        t[0]=-1.0;
        mrep['T']=t;
        t[0]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(6)
        t.clear();
        t.push_back(-1.5);
        mrep['A']=t;
        t[0]=-0.5;
        mrep['C']=t;
        t[0]=0.5;
        mrep['G']=t;
        t[0]=1.5;
        mrep['T']=t;
        t[0]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(7)
        t.clear();
        t.push_back(0.0);
        t.push_back(0.0);
        t.push_back(1.0);
        mrep['A']=t;
        t[0]=-2.0*sqrt(2.0)/3.0;
        t[1]=sqrt(6.0)/3.0;
        t[2]= -1.0/3.0;
        mrep['C']=t;
        t[0]=-2.0*sqrt(2.0)/3.0;
        t[1]=-sqrt(6.0)/3.0;
        t[2]= -1.0/3.0;
        mrep['G']=t;
        t[0]=2.0*sqrt(2.0)/3.0;
        t[1]=0;
        t[2]= -1.0/3.0;
        mrep['T']=t;
        t[0]=0;
        t[1]=0;
        t[2]=0;
        mrep['N']=t;
        mreps.push_back(mrep);


        //(8)
        t.clear();
        t.push_back(0.0);
        t.push_back(0.0);
        t.push_back(1.0);
        t.push_back(0.0);
        mrep['A']=t;
        t[0]=1.0;
        t[1]=0;
        t[2]= 0.0;
        t[3]=0.0;
        mrep['C']=t;
        t[0]=0;
        t[1]=1.0;
        t[2]= 0.0;
        t[3]=0.0;
        mrep['G']=t;
        t[0]=0;
        t[1]=0;
        t[2]= 0.0;
        t[3]=1.0;
        mrep['T']=t;
        t[0]=0;
        t[1]=0;
        t[2]= 0.0;
        t[3]=0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        //(9)
        t.clear();
        t.push_back(1.0);
        t.push_back(-1.0);
        t.push_back(-1.0);
        mrep['A']=t;
        t[0]=1.0;
        t[1]=1.0;
        t[2]=1.0;
        mrep['C']=t;
        t[0]=-1.0;
        t[1]=1.0;
        t[2]= -1.0;
        mrep['G']=t;
        t[0]=-1.0;
        t[1]=-1.0;
        t[2]= 1.0;
        mrep['T']=t;
        t[0]=0.0;
        t[1]=0.0;
        t[2]= 0.0;
        mrep['N']=t;
        mreps.push_back(mrep);
        
        mrep.clear();
        t.clear();
    }
    
    //create reps for all the contigs
    void create_reps(std::vector <std::string> nseq, std::vector <std::vector <float>> *feat,std::vector <int> *len, int dimy,std::string method)
    {
        len->push_back((int)nseq[0].size());
        std::vector <float> tmp;
        std::vector <float>  tt;
        for (int i=0;i<nseq.size();i++)
        {
            std::string t=nseq[i];
            if((int)nseq[i].size()<len->at(0))
                len->at(0)=(int)nseq[i].size();
            for (int j=0;j<nseq[i].size();j++)
            {
                for (int ii=0;ii<dimy;ii++)
                {
                    tmp.push_back(0);
                }
            }
            for (int j=0;j<nseq[i].size();j++)
            {
                int ind=repsl[method];
                auto search = mreps[ind].find(t[j]);
                tt=search->second;
                for (int ii=0;ii<dimy;ii++)
                {
                    tmp[ii*nseq[i].size()+j]=tt[ii];
                }
                tt.clear();
                std::vector<float>().swap(tt);
            }
            feat->push_back(tmp);
            tmp.clear();
            std::vector<float>().swap(tmp);
        }
    }
    
    //reps for just one contige **this is used in the main
    void create_reps(std::string nseq, std::vector <float> *feat,int *len, int dimy,std::string method)
    {
            len[0]=(int)nseq.size();
            for (int j=0;j<nseq.size();j++)
            {
                for (int ii=0;ii<dimy;ii++)
                {
                    feat->push_back(0);
                }
            }
        
            std::vector <float> ttt;
            for (int j=0;j<nseq.size();j++)
            {
                int ind=repsl[method];
                auto search = mreps[ind].find(toupper(nseq[j]));
                if((search->second).size()>0 && (search->second).size()<=10)
                {
                ttt=search->second;
                    for (int ii=0;ii<dimy;ii++)
                    {
                        feat->at(ii*nseq.size()+j)=ttt[ii];
                    }
                }
                ttt.clear();
                std::vector<float>().swap(ttt);
            }
    }

    ~nucreps()
    {
        mreps.clear();
        repsl.clear();
    }
};
#endif /* nucreps_h */
