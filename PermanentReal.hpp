//
//  PermanentReal.hpp
//  GBSNoncollision
//
//  Created by WuBujiao on 8/1/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#ifndef PermanentReal_hpp
#define PermanentReal_hpp

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>

#define MAXY 120

int IntPow(int a, int b);
int SetValueofX(int num, int x[], int t);
int GrayChange(int n, int s);
int BinaryToGray(int s);
int Permanent(int A[][MAXY], int n);

#endif /* PermanentReal_hpp */
