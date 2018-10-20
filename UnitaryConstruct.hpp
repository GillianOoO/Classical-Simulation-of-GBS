//
//  UnitaryConstruct.hpp
//  Gurvits
//
//  Created by WuBujiao on 3/12/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#ifndef UnitaryConstruct_hpp
#define UnitaryConstruct_hpp


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

/*the threshold of matrix A.*/
#define MAXZ 120

struct Norm{
    double normR;
    double normI;
};

/* 随机的产生一个n*2n的矩阵来表示n*n的复矩阵 */
void randMatrix(int n, double A[][MAXZ]);

/* 返回p和q的内积 */
Norm Inner(double p[], double q[], int n);

/*产生一个n*2n的矩阵A，其代表的n*n矩阵是一个酉矩阵（复数）*/
void SchmidtOrth(int n, double A[][MAXZ]);
#endif /* UnitaryConstruct_hpp */
