//
//  kmeans++.h
//  cluster_matagenom_nt2aa
//
//  Created by Samaneh Kouchaki on 11/07/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef kmeans___h
#define kmeans___h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
typedef struct { double x, y; int group; } point_t, *point;
#	define for_len for (j = 0, p = pts; j < len; j++, p++)
#	define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
class kmeanspp
{
public:
double randf(double m)
{
    return m * rand()/ (RAND_MAX - 1.) ;
}
    
point gen_xy(double *data ,int dim)
{
    point p, pt =(point) malloc(sizeof(point_t) * dim);
    
    int i=0;
    for (p = pt + dim; p-- > pt;) {
        p->x = data[2*i];
        p->y = data[2*i+1];
        
        i++;
    }
    return pt;
}

double dist2(point a, point b)
{
    double x = a->x - b->x, y = a->y - b->y;
    return sqrt(x*x + y*y);
}

int nearest(point pt, point cent, int n_cluster, double *d2)
{
    int i, min_i=0;
    point c;
    double d, min_d=10000000;
    

    for_n {
        min_d = 100000;
        min_i = pt->group;
        for_n {
            d = dist2(c, pt);
         //   std::cout<<d<<std::endl;
            if (min_d > (d)) {
                min_d = d; min_i = i;
            }
        }
    }
    if (d2) *d2 = min_d;
    return min_i;
}

void kpp(point pts, int len, point cent, int n_cent)
{

    int  j;
    int n_cluster;
    double sum, *d = (double *)malloc(sizeof(double) * len);
    
    point p;
    cent[0] = pts[ rand() % len ];
    for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
        sum = 0;
        for_len {
            nearest(p, cent, n_cluster, d + j);
            sum += d[j];
        }
        sum = randf(sum);
        for_len {
            if ((sum -= d[j]) > 0) continue;
            cent[n_cluster] = pts[j];
            break;
        }
    }
    for_len p->group = nearest(p, cent, n_cluster, 0);
    free(d);
}

    //clustering using kmeans++
point kmeans(point pts, int len, int n_cluster)
{
    int i, j, min_i;
    int changed;
    
    point cent = (point)malloc(sizeof(point_t) * n_cluster), p, c;
    
    //inialisation
    kpp(pts, len, cent, n_cluster);
    
    do {
        /* group element for centroids are used as counters */
        
        for_n { c->group = 0; c->x = c->y = 0; }
        for_len {
            c = cent + p->group;
            c->group++;
            c->x += p->x; c->y += p->y;
        }
        for_n { c->x /= c->group; c->y /= c->group; }
        changed = 0;
        /* find closest centroid of each point */
        for_len {
            min_i = nearest(p, cent, n_cluster, 0);
            if (min_i != p->group) {
                changed++;
                p->group = min_i;
            }
        }
    } while (changed > 0); /* stop when 99.9% of points are good */
    
    for_n { c->group = i; }
    return cent;
}

    
    //save labels and separate fasta for each cluster
void save_labels(point pts, int len, std::string name, std::vector <int> *ind, std::vector <std::string> ndata, std::vector <std::string> uids, int nparts)
{
    std::vector <std::vector <int>> final;
    std::vector <int> tmp;
    point p=pts;
    int min_i;
    int j;
    
    for (j=0;j<nparts+1;j++)
    {
        final.push_back(tmp);
    }
    
    for_len
    {
        ind->push_back(p->group+1);
    }
    std::reverse(ind->begin(),ind->end());

    std::ofstream myfile (name);
    for (j=0;j<ind->size();j++)
    {
        final[ind->at(j)].push_back(j);
        myfile<<ind->at(j)<<", ";
    }
    myfile.close();
    
    for (j=1;j<final.size();j++)
    {
        tmp.clear();
        copy((final[j]).begin(),final[j].end(),back_inserter(tmp));
        stringstream ss;
        ss << j;
        string str = ss.str();
        std::string name1=name+str;//static_cast<ostringstream*>( &(ostringstream() << j) )->str();
        std::ofstream myfile (name1);
        for (int k=0;k<tmp.size();k++)
        {
            myfile<<'>'<< uids[tmp[k]]<<"\n";
            myfile<< ndata[tmp[k]]<<"\n";
        }
        myfile.close();
    }
    final.clear();
    tmp.clear();
}
};
#endif /* kmeans___h */
