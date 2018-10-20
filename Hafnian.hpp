//
//  Hafnian.hpp
//  GBSNoncollision
//
//  Created by WuBujiao on 6/15/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#ifndef Hafnian_hpp
#define Hafnian_hpp

#include <iostream>
#include <cstdio>
#include "CHafnian.hpp"

struct Complex{
    double real;
    double image;
};

#define MAXX 120


Complex MultiComplex(Complex za, Complex zb);
Complex Hafnian(double W[][MAXX], int x[], int sa[], int n);
#endif /* Hafnian_hpp */
