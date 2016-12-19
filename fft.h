//
//  fft.h
//  cluster_metagenom
//
//  Created by Samaneh Kouchaki on 25/04/2016.
//  Copyright © 2016 Samaneh Kouchaki. All rights reserved.
//

#ifndef fft_h
#define fft_h

#include <math.h>
#include <vector>
#define PI	M_PI
#define TWOPI	(2.0*PI)

class fft
{
public:
    
    //calculate fft not used in this version
    void FFT(std::vector <float> *data, int nn, int isign)
    {
        int n, mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta;
        double tempr, tempi;
        
        n = nn << 1;
        j = 1;
        for (i = 1; i < n; i += 2) {
            if (j > i) {
                tempr = data->at(j);     data->at(j) = data->at(i);     data->at(i) = tempr;
                tempr = data->at(j+1); data->at(j+1) = data->at(i+1); data->at(i+1) = tempr;
            }
            m = n >> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }
        mmax = 2;
        while (n > mmax) {
            istep = 2*mmax;
            theta = TWOPI/(isign*mmax);
            wtemp = sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (m = 1; m < mmax; m += 2) {
                for (i = m; i <= n; i += istep) {
                    j =i + mmax;
                    tempr = wr*data->at(j)   - wi*data->at(j+1);
                    tempi = wr*data->at(j+1) + wi*data->at(j);
                    data->at(j)   = data->at(i)   - tempr;
                    data->at(j+1) = data->at(i+1) - tempi;
                    data->at(i) += tempr;
                    data->at(i+1) += tempi;
                }
                wr = (wtemp = wr)*wpr - wi*wpi + wr;
                wi = wi*wpr + wtemp*wpi + wi;
            }
            mmax = istep;
        }
    }};

#endif /* fft_h */
