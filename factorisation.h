//
//  factorisation.h
//  cluster_matagenom_nt2aa
// rsvd
//  Created by Samaneh Kouchaki on 02/06/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef factorisation_h
#define factorisation_h

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

using namespace Eigen;

class factorise
{
    MatrixXf data;
    MatrixXf transformed;
    
public:
    
    //make a matrix from the data to be factorised
    void unfold(std::vector <std::vector <float>> *d)
    {
        data.resize(d->size(),(d->at(0)).size());
        std::vector <float>tfold;
        
        for (int ii=0;ii<d->size();ii++)
        {
            copy((d->at(ii)).begin(),(d->at(ii)).end(),back_inserter(tfold));
            for (int jj=0;jj<tfold.size();jj++)
                data(ii,jj) = tfold[jj];
            tfold.clear();
        }
        d->clear();
    }
    
    void randsvd(int order)//randomised svd to perform svd faster
    {
        int p=std::min(order,(int)data.cols());
        MatrixXf y;
        y.resize((int)data.rows(), p);
        
        for (int i=0;i<p;i++)
        {
            MatrixXf rand=MatrixXf::Random((int)data.cols(),1);
            MatrixXf y1=data * rand;
            y.block(0,i,(int)data.rows(),1)=y1;
        }
        
        HouseholderQR<MatrixXf> qry(y);
        MatrixXf thinQ(MatrixXf::Identity(y.rows(),p));
        
        y.resize(0,0);
        
        MatrixXf T= (qry.householderQ()* thinQ).transpose() * data ;
    
        thinQ.resize(0, 0);
        
        JacobiSVD<MatrixXf> svd(T,  ComputeThinV);//calculate svd
        
        T.resize(0,0);
        
        VectorXf eigenvalues = svd.singularValues().real();
        MatrixXf V = svd.matrixV().real();
        
        std::vector<std::pair<double,VectorXf> > eigen_pairs;
        for(unsigned int i = 0; i < V.cols(); i++){
            double norm = V.col(i).norm();
            V.col(i) /= norm;
            eigen_pairs.push_back(std::make_pair(eigenvalues(i),V.col(i)));
        }
        
        sort(eigen_pairs.begin(),eigen_pairs.end(), [](const std::pair<double,VectorXf> a, const std::pair<double,VectorXf> b) -> bool {return (a.first > b.first);} );
        V.conservativeResize(V.rows(),p);
        
        for(unsigned int i = 0; i < p; i++)
        {
            V.col(i) = eigen_pairs[i].second;
        }
        
        eigen_pairs.clear();
        transformed = data * V;
        data.resize(0,0);
        V.resize(0,0);
    }

    //return final data after data reduction
    void get_transformed(std::vector <std::vector <float>> *result)
    {
        std::vector <float> tmp;
        for (int ii=0;ii<transformed.rows();ii++)
        {
            for (int jj=0;jj<transformed.cols();jj++)
            {
                tmp.push_back(transformed(ii,jj));
            }
            result->push_back(tmp);
            tmp.clear();
        }
        std::vector<float>().swap(tmp);
        transformed.resize(0,0);
    }
    
};
#endif /* factorisation_h */
