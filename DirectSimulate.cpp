//
//  DirectSimulate.cpp
//  GBS
//
//  Created by WuBujiao on 5/30/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#include "DirectSimulate.hpp"


/* 随机生成一个 0~1 之间的随机数 */
double Myrand(){
    
    double r;
    //srand(unsigned(time(NULL)));
    r = rand()*1.0/RAND_MAX;
    return r;
}
/* Generate a sample of given probability distribution (p_0,...,p_{n-1}), 返回N,表示以概率p(N)得到 N */
int DirectSimulate(double p[]){
    int N=0;
    //double p[100];
    double F=0;
    double U;
    U=Myrand();
  //  std::cout<<"U"<<U<<' ';
    while(F<U) {
        F += p[N++];
    }
    N--;
  //  std::cout<<N<<std::endl;
    return N;
}
