//
//  Hafnian.cpp
//  GBSNoncollision
//
//  Created by WuBujiao on 6/15/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#include "Hafnian.hpp"


Complex MultiComplex(Complex za, Complex zb){
    Complex z;
    z.real = za.real * zb.real - za.image * zb.image;
    z.image = za.real * zb.image + za.image * zb.real;
    return z;
}

Complex Hafnian(double W[][MAXX], int x[], int s[], int l){
    Complex z;
    if (l == 0) {
        z.image = 0;
        z.real = 1;
    }
    if (l == 2) {
        z.real = W[x[s[0]]][x[s[1]]*2];
        z.image = W[x[s[0]]][x[s[1]]*2+1];
        
    }//W_{1,2}
    
    Complex za,zb;
    if (l == 4) {
        int i,j;
        i = 0,j=1;
        za.real = W[x[s[i]]][x[s[j]]*2];
        za.image = W[x[s[i]]][x[s[j]]*2+1];
        i = 2, j = 3;
        zb.real = W[x[s[i]]][x[s[j]]*2];
        zb.image = W[x[s[i]]][x[s[j]]*2+1];
        z = MultiComplex(za, zb);
        
        i = 0,j=2;
        za.real = W[x[s[i]]][x[s[j]]*2];
        za.image = W[x[s[i]]][x[s[j]]*2+1];
        i = 1, j = 3;
        zb.real = W[x[s[i]]][x[s[j]]*2];
        zb.image = W[x[s[i]]][x[s[j]]*2+1];
        z.real += MultiComplex(za, zb).real;
        z.image += MultiComplex(za, zb).image;
        
        i = 0,j=3;
        za.real = W[x[s[i]]][x[s[j]]*2];
        za.image = W[x[s[i]]][x[s[j]]*2+1];
        i = 1, j = 2;
        zb.real = W[x[s[i]]][x[s[j]]*2];
        zb.image = W[x[s[i]]][x[s[j]]*2+1];
        z.real += MultiComplex(za, zb).real;
        z.image += MultiComplex(za, zb).image;
    }
    
    if (l >= 6) {
        double complex a[l*l];
        int lj = 0;
        for (int j = 0; j < l; j ++) {
            for (int k = 0; k < l; k ++) {
                //a[j][k]={ W[x[s[j]]][x[s[k]]*2], W[x[s[j]]][x[s[k]] * 2 + 1]};
                a[lj ++] = W[x[s[j]]][x[s[k]]*2] + I * W[x[s[j]]][x[s[k]] * 2 + 1];
            }
        }
        double complex fz;
        /*
         call hafnian library in this place,
         input:
         a: a l*l even square matrix.
         l: the size of a.
         Return : fz = Hafnian(a,l);
         */
        fz = lhafnian(a, l);
        z.real = creal(fz);
        z.image = cimag(fz);
    }
    return z;

}