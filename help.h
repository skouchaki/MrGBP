//
//  help.h
//  cluster_visualise_stft
//
//  Created by Samaneh Kouchaki on 04/08/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef help_h
#define help_h

#include <stdio.h>

class help
{
public:
    help()
    {
        std::cout<<"To run the code you need to indicate the input file. There are several other optional parameters that can be changed. All options and their default values are as follows: " <<std::endl;
        std::cout<<"-fa <filepath> : indicate the fasta file to be processed by MLBP_BIN"<<std::endl;
        std::cout<<"-svdd <dim> : indicate the SVD dimension, default: 60"<<std::endl;
        std::cout<<"--mlbpn <dim> : indicate the MLBP window length, default: 8"<<std::endl;
        std::cout<<"-outdir <filepath> : indicate a directory that all the MLBP_BIN results will be saved including cluster labels and .png figure"<<std::endl;
        std::cout<<"-no_clust <num> : indicate the number of clusters for kmeans++, default: 10"<<std::endl;
        std::cout<<"-covpm <num> : indicate a file contains average coverage depth"<<std::endl;
        std::cout<<"-covpstd <num> : indicate a file contains standard deviation coverage depth"<<std::endl;
        std::cout<<"-reps <method> : indicate the the numerical mapping including EIIP, Real, Integer, Paired, Atomic, default: EIIP"<<std::endl;
        std::cout<<"-clust <method> : indicate the clustering method kmeans++ or dbscan: dbscan"<<std::endl;
        std::cout<<"-dbep <num> : indicate the epsilon for dbscan, default: 0.02"<<std::endl;
        std::cout<<"-dbminpt <num> : indicate the number of minimum neighbouring points for dbscan, default: 8"<<std::endl;
        std::cout<<"-save_feat : If the features needs to be saved separately"<<std::endl;
        std::cout<<"-no_dims : bh-tSNE dimensions, defualt: 2"<<std::endl;
    }
    
};
#endif /* help_h */
