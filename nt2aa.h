//
//  nt2aa.h
//  cluster_matagenom_nt2aa
//
//  Created by Samaneh Kouchaki on 05/07/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef nt2aa_h
#define nt2aa_h

#include <string>
#include <iostream>
#include <vector>

using namespace std;

class nt2aa
{
public:
    //nucl to amino acid not used in this version
    void con_fasta(std::vector <std::string> ndata, std::vector <std::string> ids, std::string name)
{
    std::ofstream myfile (name.c_str());
    if (myfile.is_open())
    {
        for (int i=0;i<ndata.size();i++)
        {
            myfile<< '>'<<ids[i]<<"\n";
            myfile<<con_nt2aa(ndata[i]);
            myfile<<"\n";
        }
        myfile.close();
    }
    else printf("Unable to open file");
}
string con_nt2aa(string ndata)
{
    string tmp1="";
    for (int j=0;j<floor(ndata.size()/3);j++)
    {
        tmp1+=convert(ndata.substr(j*3,3));
    }
 //   cout <<tmp1<<endl;
    return tmp1;
}

char convert(string input)
{
    char aa='*';
    switch(input[0])
    {
        case 'A':
        switch(input[1])
        {
            case 'A':
            switch(input[2])
            {
                case 'A':
                case 'G':
                    aa='K';
                    break;
                case 'C':
                case 'T':
                    aa='N';
                    break;
            }
            break;
            case 'C':
            switch(input[2])
            {
                case 'A':
                case 'G':
                case 'C':
                case 'T':
                    aa='T';
                    break;
            }
            break;
            case 'G':
            switch(input[2])
            {
                case 'A':
                case 'G':
                    aa='R';
                    break;
                case 'C':
                case 'T':
                    aa='S';
                    break;
            }
            break;
            case 'T':
            switch(input[2])
            {
                case 'G':
                    aa='M';
                    break;
                case 'A':
                case 'C':
                case 'T':
                    aa='I';
                    break;
                    }
            }
            break;
            
        
            case 'C':
            switch(input[1])
            {
                case 'A':
                    ;
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                        aa='Q';
                        break;
                    case 'C':
                    case 'T':
                        aa='H';
                        break;
                }
                break;
                case 'C':
                switch(input[2])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        aa='P';
                        break;
                }
                break;
                case 'G':
                switch(input[2])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        aa='R';
                        break;
                }
                break;
                case 'T':
                switch(input[2])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        aa='L';
                        break;
                }
            }
            break;
            
        
            case 'G':
            switch(input[1])
            {
                case 'A':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                        aa='E';
                        break;
                    case 'C':
                    case 'T':
                        aa='D';
                        break;
                }
                break;
                case 'G':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                    case 'C':
                    case 'T':
                        aa='G';
                        break;
                }
                break;
                case 'C':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                    case 'C':
                    case 'T':
                        aa='A';
                        break;
                }
                break;
                case 'T':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                    case 'C':
                    case 'T':
                        aa='V';
                        break;
                }
                break;
            }
            break;
        
        
            case 'T':
            switch(input[1])
            {
                case 'A':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                        aa='*';
                        break;
                    case 'C':
                    case 'T':
                        aa='Y';
                        break;
                }
                break;
                case 'G':
                switch(input[2])
                {
                    case 'A':
                        aa='*';
                        break;
                    case 'G':
                        aa='W';
                        break;
                    case 'C':
                    case 'T':
                        aa='C';
                        break;
                }
                break;
                case 'C':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                    case 'C':
                    case 'T':
                        aa='S';
                        break;
                }
                break;
                case 'T':
                switch(input[2])
                {
                    case 'A':
                    case 'G':
                        aa='L';
                        break;
                    case 'C':
                    case 'T':
                        aa='F';
                        break;
                }
                break;
            }
    }
    return aa;
}
    
};
#endif /* nt2aa_h */
