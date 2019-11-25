//
//  Bruteforce.cpp
//  ExpGBS
//
//  Created by WuBujiao on 12/22/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#include "Bruteforce.hpp"


int FAC[100];

double f;

int flag;

//int count;

// z_1, ..., z_k-1 given, compute z_k and dfs ...
void dfs(int z[MAXX], int k, int n, double x, double **A){
    if(k == n || flag == 1)
    return;
//    count ++;
    int s[100], sa[100];
    double ans;
   // for(int j = 0; j < n * n; j++){
   //    s[j] = 0;
   // }
   // for (int j = 0; j < n; j++)sa[j] = j;

    for (int j = z[k-1]; j < n * n; j ++) {
        z[k] = j;
//        for(int j = 0; j < n * n; j++){
//                s[j] = 0;
//        }
//        for(int i = 0; i < n; i ++){
//                s[z[i]]++;
//        }

        dfs(z, k + 1, n, x, A);
	for(int j = 0; j < n * n; j++){
       		s[j] = 0;
    	}
	for(int i = 0; i < n; i ++){
		s[z[i]]++;
	}
//	for (int i = 0; i < n; i++)cout<<z[i]<<" ";
//	for (int i = 0; i < n * n; i++){
//		cout<<s[i]<<" ";
//	}
//	cout<<endl;
   	for (int j = 0; j < n; j++)sa[j] = j;
	Complex com = Hafnian(A, z, sa, n);
	ans = com.real * com.real + com.image * com.image;
	for (int i = 0; i < n * n; i++){
		ans /= FAC[s[i]];
//		if(FAC[s[i]] == 0)cout<<"Error! "<<s[i]<<endl;
	//	cout<<s[i]<<": "<<FAC[s[i]]<<" and "<<ans<<endl;
	}
	f += ans;
//	cout<<f<<endl;
	if(f > x && flag == 0){
		flag = 1;
		/*output the answer of sampling.*/
	//	for (int i = 0; i < n; i++){
	//		cout<<z[i]<<"\t";
	//	}
	//	cout<<endl;
	}
      //  for (int i = 0; i < n; i ++) {
      //      cout<<z[i]<<"\t";
      //  }
      //  cout<<endl;
    }
    
    
}

double B(int n, int k){
    double num = 1;
    double dec = 1;
    if(n < k)return 0;
    if(n - k < k)k = n - k;
    for(int j = 0; j < k; j++){
        num *= n-j;
        dec *= k-j;
    }
    return num/dec;
}

void bruteforce(int n, double **A){
    
    //FILE *fp;
   // fp = fopen("Bruteforce.txt", "w+");
           
    int z[MAXX];
    double x = rand() * 1.0/RAND_MAX;
    //cout<<x<<endl;
   // x = 1;
    x = x * B((n * n + n)/2 - 1, n/2);
//    cout<<x<<endl;

    f = 0;
    flag = 0;
    FAC[0] = 1;
//    count = 0;
    for (int j = 1; j <= n; j++){
	FAC[j] = FAC[j - 1]* j;
    }
    for (int j = 0; j < n * n; j++){
        z[0] = j;
	dfs(z, 1, n, x, A);

    } 
    
   // cout<<"length: "<<count<<endl;
   // fclose(fp);
    
}

