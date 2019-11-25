//
//  Bruteforce.hpp
//  ExpGBS
//
//  Created by WuBujiao on 12/22/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#ifndef Bruteforce_hpp
#define Bruteforce_hpp

#include <stdio.h>
#include <iostream>
#include "UnitaryConstruct.hpp"
#include "Hafnian.hpp"
#include "DirectSimulate.hpp"


using namespace std;
//int p[10000];

void bruteforce(int n, double **A);

void dfs(int z[MAXX], int k, int n, double x, double **A);
#endif /* Bruteforce_hpp */
