//
//  load.h
//  cluster_matagenom_nt2aa
//
//  Created by Samaneh Kouchaki on 10/05/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef load_h
#define load_h

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>

class load
{
public:
    
    //reverse_complement
    std::string reverse_com(std::string inp)
    {
        std::string rev;
        for (int i=0;i<inp.size();i++)
        {
            switch(inp[i])
            {
                case 'A':
                    rev+='T';
                    break;
                case 'T':
                    rev+='A';
                    break;
                case 'G':
                    rev+= 'C';
                    break;
                case 'C':
                    rev+= 'G';
                    break;
                case 'N':
                    rev+= 'N';
                    break;
            }
        }
        return rev;
    }
    
    
    //load a fasta file indicated by the directory
    load(std::string filepathn,std::vector <std::string> * ndata, std::vector <std::string> * uids)
    {
        ////////nucloetide
        std::string uid, read;
        int count = 0;
        
        std::string lineentry, tmp;
        std::ifstream ngsFile1(filepathn.c_str());
        
        int i=0;
        int j=0;
        //Start parsing the fasta  file until it reach the end of the file
        while(!ngsFile1.eof())
        {
            getline(ngsFile1, lineentry);
            std::stringstream readstream(lineentry);
            
            if(count == 0 && lineentry[0] == '>')
            {
                //If line starts with > char and count = 0 then asinge line content to uid
                //1. save previous read
                if(i>0 )
                {
                    std::size_t found = uid.find("complement");
                    if (found==std::string::npos)
                    {
                        ndata->push_back(read);
                        uids->push_back(uid);
                    }
                    else
                    {
                        read=reverse_com(read);
                        read.erase(std::remove(read.begin(), read.end(), 'N'), read.end());
                        if(read.size()>=1000)
                        {
                        ndata->push_back(read);
                        uids->push_back(uid);
                        }
                    }
                }
                i++;
                //2. new id
                uid = lineentry;
                uid.erase(0,1);
                count++;
                read = "";
            }
            else if (count ==0)//remaining of the current read
            {
                tmp  = "";
                readstream >> tmp;
                read=read+tmp;
            }
            else
            {
                //If line doesn't starts with > char and count = 1  asinge line content to read
                readstream >> read;
                count=0;
            }
        }
        //Close file
        ngsFile1.close();
    }
};
#endif /* load_h */
