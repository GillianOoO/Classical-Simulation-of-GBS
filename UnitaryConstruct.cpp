//
//  UnitaryConstruct.cpp
//  Gurvits
//
//  Created by WuBujiao on 3/12/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#include "UnitaryConstruct.hpp"



double MyRand(){
    return rand()*1.0/RAND_MAX;
}
/* 随机的产生一个n*2n的矩阵来表示n*n的复矩阵 */
void randMatrix(int n, double A[][MAXZ]){
    int i,j;
    
    for (i = 0; i<n; i++) {
        for (j = 0; j<2*n; j++) {
            A[i][j] = MyRand();
        }
    }
}

/* 返回p和q的内积 */
Norm Inner(double p[], double q[], int n){
    // double normR=0,normI = 0;
    Norm norm;
    norm.normR = 0;
    norm.normI = 0;
    for (int i=0; i<n; i++) {
        norm.normR += p[2*i]*q[2*i] + p[2*i+1]*q[2*i+1];
        norm.normI += p[2*i]*q[2*i+1] - p[2*i+1]*q[2*i];
    }
    return norm;
    
}

/*产生一个n*2n的矩阵A，其代表的n*n矩阵是一个酉矩阵（复数）*/
void SchmidtOrth(int n, double A[][MAXZ]){
    
    Norm norm,inn;
    double z[MAXZ];
    srand((unsigned)time(NULL));
    randMatrix(n, A);
    norm = Inner(A[0],A[0],n);
    
    for (int i=0; i<2*n; i++) {
        A[0][i] /= sqrt(norm.normR);
    }
    
    for (int j = 1; j < n; j++ ) {
        for (int i = 0; i< 2*n; i++) {
            z[i] = 0;
        }
        //        memset(z, 0, 2*n);
        
        for (int i = 0; i < j; i++) {
            inn = Inner(A[i],A[j],n);  //Inner = <A_i |A_j >
            for (int k = 0; k < n; k++) {
                z[2*k] += A[i][2*k]*inn.normR - A[i][2*k+1]*inn.normI;
                z[2*k+1] += A[i][2*k]*inn.normI + A[i][2*k+1]*inn.normR;
            }// Inner * |A_i>
        }//sum_{j<i}<A_i|A_j>*|A_i>
        
        for (int k = 0; k < 2*n ; k++) {
            A[j][k] -= z[k];
        }
        norm = Inner(A[j],A[j],n);
        for (int k = 0; k < 2*n ; k++) {
            A[j][k] /= sqrt( norm.normR);
        }
    }
}