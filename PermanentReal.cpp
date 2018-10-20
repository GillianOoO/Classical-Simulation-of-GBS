//
//  PermanentReal.cpp
//  GBSNoncollision
//
//  Created by WuBujiao on 8/1/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#include "PermanentReal.hpp"

int IntPow(int a, int b){
    int k, power=1;
    for (k=0; k<b; k++) {
        power*=a;
    }
    return power;
}
/* 将num转化为二进制的t位数，存在x数组中，并返回奇偶性:偶数：1， 奇数：-1。 */
int SetValueofX(int num, int x[], int t){
    int parity, flag=0;
    for (int i=0; i<t; i++) {
        x[i]= num%2;
        flag += x[i];
        num /= 2;
        //    cout<<x[i];
    }
    //  cout<<endl;
    parity = (flag%2 == 0? 1:-1);//(-1)^(x_1 + ... + x_n)
    return parity;
}
/* 返回s和s-1对应glay码的不相同的比特位k,若是由0 --> 1: return -1
 若是由1 --> 0: return 1。
 */
int GrayChange(int n, int s){
    int k=0;
    
    s--;
    while (s%2==1) {
        s /= 2;
        k++;
    }
    //判断格雷码g-1 --> g是由0-->1还是由1-->0.若是由0-->1返回-k，否则返回+k
    s /= 2;
    k = (k== 0? n: k);
    return s % 2==0? -1*k : k;
}

int BinaryToGray(int s){
    return s ^ (s>>1);
}

/*A为n*n实矩阵。*/
int Permanent(int A[][MAXY], int n){
    int x[MAXY];
    int factor, flag, fflag;
    double sum[MAXY], perm, temp=0;
    //memset(sum, 0, n);
    
    for (int i = 0; i < n; i ++) {
        sum[i] = 0;
    }
    if (n == 0) {
        return 1;
    }
    for (int j = 0; j< n; j++) {
        // sum = 0;
        for (int k = 0; k < n; k++) {
            sum[j] += A[j][k];
        }
        temp = (j == 0? sum[j] : temp * sum[j]);
    }
    perm = temp;
    for (int i=1; i < IntPow(2, n-1); i++) {
        factor = SetValueofX(BinaryToGray(i),x,n);
        flag = GrayChange(n, i);
        fflag = abs(flag) == n ? 0 : abs(flag);
        for(int j=0;j<n;j++){
            sum[j] += flag/abs(flag) * 2 * A[j][fflag];
            temp = (j==0? sum[j] : (temp * sum[j]));
            //        cout<<j<<" "<<temp<<endl;
        }
        //  cout<<endl;
        perm += temp * factor;
        
    }
    perm /= IntPow(2, n-1);
    
    return perm;
    
}