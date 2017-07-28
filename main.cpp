//
//  main.cpp
//
//  Created by Samaneh Kouchaki on 03/08/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

//#pragma GCC diagnostic ignored "-Wreturn-type-c-linkage"
// g++ *.cpp -Wno-return-type-c-linkage -I /usr/local/Cellar/eigen/3.2.8/include/eigen3/  -o a  -I /usr/local/include/ -lmgl -std=gnu++11
//./a -fa /Users/samaneh/Downloads/bowtie2-2.2.9/refre/sorted_contigs.fa -covpm /Users/samaneh/Downloads/bowtie2-2.2.9/refre/coverage_profile  -covpstd /Users/samaneh/Downloads/bowtie2-2.2.9/refre/Coverage_std -outdir /Users/samaneh/Downloads/bowtie2-2.2.9/refre/sortedres

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <map>
#include <sys/stat.h>

#include "load.h"
#include "help.h"
#include "nucreps.h"
#include "feature_extraction.h"
#include "factorisation.h"
#include "tsne.h"
#include "kmeanspp.h"
#include "dbscan.h"

#include <mgl2/mgl.h>

//saving BH-tSNe as .png file
int sample(mglGraph *gr,double *a,int dim,std::vector <int> id,int k,int no_dims)
{
    
    mglData x;
    mglData y;
    mglData z;
    double * y1=new double [dim];
    double min=100000;
    double max=-100000;
    for (int i=0;i<dim;i++)
    {
        y1[i]=a[no_dims*i];
        if(y1[i]>max)
            max=y1[i];
        else if(y1[i]<=min)
            min=y1[i];
    }
    
    x.Set(y1,dim);
    gr->SetRange('x', min-1, max+1);
    
    min=100000;
    max=-100000;
    
    for (int i=0;i<dim;i++)
    {
        y1[i]=a[no_dims*i+1];
        if(y1[i]>max)
            max=y1[i];
        else if(y1[i]<=min)
            min=y1[i];
    }
    y.Set(y1,dim);
    gr->SetRange('y', min-1, max+1);
    
    if(no_dims==3)
    {
        min=100000;
        max=-100000;
        
        for (int i=0;i<dim;i++)
        {
            y1[i]=a[no_dims*i+1];
            if(y1[i]>max)
                max=y1[i];
            else if(y1[i]<=min)
                min=y1[i];
        }
        z.Set(y1,dim);
        gr->SetRange('z', min-1, max+1);
    }
    gr->SetOrigin(0,0,0);
    gr->SetFontSize(2);
    gr->SubPlot(2,1,0,"");  gr->Title("Visualisation"); if(no_dims==3)gr->Rotate(60,40);
    gr->Box();  if(no_dims==3) gr->Plot(x,y,z," .");else  gr->Plot(x,y," .");
    
    //color printing
    gr->SubPlot(2,1,1,"");  gr->Title("Cluster Labels  Visualisation");if(no_dims==3)gr->Rotate(60,40);
    for (int j=0;j<=k;j++)
    {
        int s=0;
        for (int i=0;i<dim;i++)
            if(id[i]==j)
            {
                y1[s]=a[no_dims*i];
                s=s+1;
            }
        x.Set(y1,s);
    
        s=0;
        for (int i=0;i<dim;i++)
            if(id[i]==j)
            {
                y1[s]=a[no_dims*i+1];
                s=s+1;
            }
        y.Set(y1,s);

        if(no_dims==3)
        {
            s=0;
            for (int i=0;i<dim;i++)
                if(id[i]==j)
                {
                    y1[s]=a[no_dims*i+1];
                    s=s+1;
                }
            z.Set(y1,s);

        }
        
        if(j==0)
             gr->Box();
        if(no_dims==3)
         gr->Plot(x,y,z," .");
        else
            gr->Plot(x,y," .");
    }
    return 0;
}

//used to find the coverage informtion of each file if the file provided
string find_infile(std::string file, std::string name,std::vector <std::string> *record)
{
    std::string s;
    int flag=0;
    std::ifstream myfile (file.c_str());
    std::vector <std::string> record1;
//std::cout<<name<<std::endl;
     while(getline( myfile, s ))
     {
         record1.clear();
         std::istringstream ss( s );
         while (ss)
         {
             std::string s;
             if (!getline( ss, s, ' ' )) break;
             record1.push_back( s );
         }
//std::cout<<"aa "<<record1[0]<<std::endl;
         if (name.compare(record1[0])==0)
         {
             flag=1;
             break;
         }
         else if (name.find(record1[0])!=std::string::npos && name[record1[0].size()]==' ')
         {
             flag=1;
             break;
         }
	else if (name.find(record1[0])!=std::string::npos && name[record1[0].size()]=='\r')
         {
             flag=1;
             break;
         }
     }
    if(flag==0)
    {
        for (int i=1;i<record1.size();i++)
        {
            record->push_back("0");
        }
	myfile.close();
        return "0";
    }
    else
    {
      //  std::cout<<record1[1]<<std::endl;
        for (int i=1;i<record1.size();i++)
        {
            record->push_back(record1[i]);
        }
	myfile.close();
        return record1[1];
    }
}

int main(int argc, const char * argv[])
{
    int k=1;
  //  std::map<std::string, int> ids;
 //   std::vector <int> nnlen;
    
    std::string name;
    std::string cpname,cpname1;
    std::string out="tsne_results1";
    std::string dire;
    
    std::string femethod = "n_LBP";
    std::string repsmethod = "Integer";
    int repsdim=1;
    int lev=6;
    int dd=60;
    int nparts=10;
    float dbep=0.02;
    int dbminpt=8;
    std::vector <std::string> ndata;
    std::vector <std::string> uids;
    
    int flag=0;
    int flag1=0;
    int flag2=0;
    int flag3=0;
    int flag4=1;
    int flag5=0;
    int maxs=1000;
    
    int no_dims=2;
    int perplexity=40;
    double theta=0.5;
    
    //loop to get information from cmd
    int iii=1;
    while (iii<argc)
    {
        if ((std::string)(argv[iii]) == "-fa")//fasta file
        {
            name = (std::string) (argv[iii + 1]);
            if (FILE *file = fopen(name.c_str(), "r"))
            {
                fclose(file);
                flag=1;
            }
            else
            {
                help h;
                break;
            }
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-help")//save features separately
        {
            help h;
            iii=iii+1;
        }
        else if ((std::string)(argv[iii]) == "-svdd")//svd dim
        {
            dd=atoi(argv[iii + 1]);
            iii=iii+2;
        }
        
        else if ((std::string)(argv[iii]) == "-mlbpn")//mlbp dim
        {
            lev=atoi(argv[iii + 1]);
            iii=iii+2;
        }

        else if ((std::string)(argv[iii]) == "-outdir")//outdir
        {
            dire = (std::string) (argv[iii + 1]);
            mkdir(dire.c_str(), 0777);
            flag1=1;
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-no_clust")//no_clust in kmeans++
        {
            nparts =atoi(argv[iii + 1]);
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-dbep")//epsilon for dbscan
        {
            dbep =atof(argv[iii + 1]);
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-dbminpt")//minpt for dbscan
        {
            dbminpt =atoi(argv[iii + 1]);
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-covpm")//average depth
        {
             cpname = (std::string) (argv[iii + 1]);
             flag2=1;
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-covpstd")//std depth
        {
            cpname1 = (std::string) (argv[iii + 1]);
            flag3=1;
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-reps")//numerical mapping
        {
            repsmethod = (std::string) (argv[iii + 1]);
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-clust")//std depth
        {
            cpname1 = (std::string) (argv[iii + 1]);
            if(cpname1.compare("kmeans++"))
                flag4=0;
            iii=iii+2;
        }
        else if ((std::string)(argv[iii]) == "-save_feat")//save features separately
        {
            flag5=1;
            iii=iii+1;
        }
        else if ((std::string)(argv[iii]) == "-no_dims")//num of dimensions (BH-tsne)
        {
            no_dims=atoi(argv[iii + 1]);
            iii=iii+2;
        }
	else if ((std::string)(argv[iii]) == "-mincl")//minimum contig length
        {
            no_dims=atoi(argv[iii + 1]);
            iii=iii+2;
        }
        
    }
    
    if(flag1==0)//if not out dir indicated it makes one
    {
        dire=name+"1";
        mkdir(dire.c_str(), 0777);
    }
    if(flag==1)
    {
        //loading data
        clock_t t1=clock();
        load l(name,&ndata, &uids);
        std::cout<<"loading times: "<<(float)(clock()-t1)/CLOCKS_PER_SEC<<std::endl;
        
        //numerical reps
        nucreps *rn=new nucreps();
        std::vector <float> repsn;
        
        //features
        feature_extraction fe;
      //  std::size_t found1;
        std::vector <std::vector <float>> features;
        
        int k = 1;
      //  std::vector <int> id;
      //  std::vector <int> lens;
        int *lenn=new int(0);
        
        //loop for extracting features
        t1=clock();
   //     std::string line;
   //     std::ifstream myfile (cpname.c_str());
   //     std::ifstream myfile1 (cpname1.c_str());
        std::vector <std::string> record;
        for (int i=0;i<ndata.size();i++) //for all the contigs
        {
            if(ndata[i].size()>=maxs)//only long contigs
            {        
	//	lens.push_back(ndata[i].size());
                ndata[i].erase(std::remove(ndata[i].begin(), ndata[i].end(), 'N'), ndata[i].end());
                record.clear();
                if(flag2==1)
                {
		//	std::cout<<cpname<<", "<<uids[i]<<std::endl;
                    find_infile(cpname, uids[i],&record);
                }
                if(flag3==1)
                {
                    find_infile(cpname1, uids[i],&record);
                }

                rn->create_reps(ndata[i], &repsn,lenn,repsdim,repsmethod);
                fe.extract(femethod, &repsn,ndata[i], lenn[0], repsdim, &features,lev,record,ndata[i].size());
                ndata[i]="";
                repsn.clear();
                
                record.clear();
	//	std::cout<<uids[i]<<std::endl;
                std::string nn = find_infile("/local/engs1758/Downloads/100s/100s.spe.txt", uids[i],&record);
     //           for (int i=0;i<record.size();i++)
     //           std::cout<<nn<<", "<<uids[i]<<std::endl;
     //               std::cout<<std::endl;
             //   std::cout<<i<<std::endl;
          //      if(ids.count(nn)==0)
         //       {
                //    std::cout<<nn<<std::endl;
           //         ids[nn]=k;
           //         k++;
           //     }
           //     auto search = ids.find(nn);
           //     id.push_back(search->second);
                
            }
        }
        std::cout<<"feature selectiom times: "<<(float)(clock()-t1)/CLOCKS_PER_SEC<<std::endl;
        
       /* std::ofstream myfi ("/Users/samaneh/Downloads/bowtie2-2.2.9/simulated_data/labs100.txt");
        if (myfi.is_open())
        {
            for (int i=0;i<id.size();i++)
            {
                myfi << id[i] <<" "<< nnlen[i] << "\n";
            }
            myfi.close();
        }
        id.clear();
        nnlen.clear();
        std::cout<<ids["0"]<<std::endl;*/
        
       // free(lenn);
        ndata.clear();
        uids.clear();
        std::vector <double> tmp;
        if(flag5)
        {
            out=dire+"/MLBP_features";
            std::ofstream myfile2 (out);
            if (myfile2.is_open())
            {
                for (int i=1;i<=features.size();i++)
                {
                    tmp.clear();
                    copy((features[i]).begin(),(features[i]).end(),back_inserter(tmp));
                    for (int j=0;j<tmp.size();j++)
                        myfile2<< tmp[j]<<" ";
                    //  myfile<<id[i]<<" "<<lens[i];
                    myfile2<<"\n";
                }
                myfile2.close();
            }
            else printf("Unable to open file");
        }
        
       // std::cout<<features.size()<<std::endl;
        //rsvd dimension redection
        t1=clock();
        std::vector<std::vector <float>> trans;
        factorise fz;
        fz.unfold(&features);
        fz.randsvd(dd);
        fz.get_transformed(&trans);
        std::cout<<"svd time: "<<(float)(clock()-t1)/CLOCKS_PER_SEC<<std::endl;
        
        //bh-tsne
        int dimx=(int)trans.size();
        int dim=dd;
        double* Y = (double*) malloc(dimx* no_dims * sizeof(double));
        double* costs = (double*) calloc(dimx, sizeof(double));
        if(Y == NULL || costs == NULL) { printf("Memory allocation failed!\n"); exit(1); }
        TSNE *tsne=new TSNE();
        
        double* data =new double[dimx*dim];
        for (int i=0;i<dimx;i++)
        {
            tmp.clear();
            copy((trans[i]).begin(),(trans[i]).end(),back_inserter(tmp));
            for (int j=0;j<tmp.size();j++)
            {
                data[i*dim+j]=tmp[j];
            }
        }
        trans.clear();
        out=dire+"/"+"tsne_results";
        
        t1=clock();
        tsne->run(data, dimx, dim, Y, no_dims, perplexity, theta,-1 , false);
	//std::cout<<"done"<<std::endl;        
	tsne->save_data(Y, dimx, no_dims,out);
        std::cout<<"dimension reduction running time: "<<(float)(clock()-t1)/CLOCKS_PER_SEC<<std::endl;
        
        free(tsne);
        free(data);
        
        //clustering
        t1=clock();
        if(flag4==0)
        {
            std::vector <int> cen;
            kmeanspp kmpp;
            point v= kmpp.gen_xy(Y,dimx);
            point ccc = kmpp.kmeans(v, dimx, nparts);
            out=dire+"/"+"labels";
            load l1(name,&ndata, &uids);
            std::vector <std::string> ndata1;
            std::vector <std::string> uids1;
            for (int i=0;i<ndata.size();i++) //for all the reads
            {
                if(ndata[i].size()>=maxs)
                {
                    ndata1.push_back(ndata[i]);
                    uids1.push_back(uids[i]);
                }
            }
            ndata.clear();
            uids.clear();
            kmpp.save_labels(v, dimx, out,&cen,ndata1,uids1,nparts);
            ndata1.clear();
            uids1.clear();
        
            //save png file
            mglGraph gr;
            gr.Alpha(true);
            gr.Light(true);
            sample(&gr,Y,dimx,cen,nparts,no_dims); // The same drawing function.
            out=dire+"/"+"figure1.png";
            gr.WritePNG(out.c_str());  // Don't forget to save the result!
        }
        else
        {
        
            clustering::DBSCAN::ClusterData cl_d = clustering::DBSCAN::gen_cluster_data( no_dims, dimx,Y );
            clustering::DBSCAN db(dbep,dbminpt,1);
        
            db.fit( cl_d );
            out=dire+"/"+"labels";
            load l1(name,&ndata, &uids);
            std::vector <int> labs;
            std::vector <std::string> ndata1;
            std::vector <std::string> uids1;
            for (int i=0;i<ndata.size();i++) //for all the reads
            {
                if(ndata[i].size()>=maxs)
                {
                    ndata1.push_back(ndata[i]);
                    uids1.push_back(uids[i]);
                }
            }
            ndata.clear();
            uids.clear();
            nparts=db.save_labels(out.c_str(),ndata1,uids1,&labs,Y);
            mglGraph gr;
            gr.Alpha(true);
            gr.Light(true);
            sample(&gr,Y,dimx,labs,nparts,no_dims); // The same drawing function.
            out=dire+"/"+"figure1.png";
            gr.WritePNG(out.c_str());  // Don't forget to save the result!
        }
        std::cout<<"clustering time: "<<(float)(clock()-t1)/CLOCKS_PER_SEC<<std::endl;
        
    }
    else
    {
        help h;
    }
    return 0;
}
