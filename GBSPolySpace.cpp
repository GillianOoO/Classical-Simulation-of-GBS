//
//  main.cpp
//  GBSNoncollision
//
//  Created by WuBujiao on 6/15/18.
//  Copyright 漏 2018 WuBujiao. All rights reserved.
//

/**
Basic classical gaussian boson sampling algorithm.
The sampling process is in GBS.


**/

#include "UnitaryConstruct.hpp"
#include "DirectSimulate.hpp"
#include "Hafnian.hpp"
#include "PermanentReal.hpp"
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <iterator>
#include <time.h>
#include <sstream>
#include "Bruteforce.hpp"

using namespace std;

#define SampleMAX 8000
#define SMAX 10000
#define MAXW 5000

double **W = NULL;

#define max(a,b) (a>b? a: b)
#define deltaZero(a) (a == 0? 0 : 1)

string to_string(int num){
    ostringstream oss;
    oss<<num;
    return oss.str();
}
/*binom{n}{k}:鑻<0鎴栬€呭ぇ浜巒杩斿洖0*/
int Binom(int n, int k){
    int ans = 1;
    // int dec = 1;
    if (k < 0 || k > n) {
        return 0;
    }
    k = (k<(n-k)? k : n-k);
    if(k == 0) return 1;
    for (int i = 1; i<= k; i++) {
        ans *= n - k + 1;
        ans /= k;
    }
    return ans;
}


/* W = AA^t W: m * 2m  */
void SymmetricM(double **A, int m){
    
    for (int i = 0; i < m ; i ++) {
        for (int j = 0; j < m; j ++) {
            for (int k = 0; k < m; k ++) {
                W[i][2*j] += A[i][2*k] * A[j][2*k] - A[i][2*k+1]*A[j][2*k+1];
                W[i][2*j+1] += A[i][2*k]*A[j][2*k+1] + A[i][2*k+1]*A[j][2*k];
            }
        }
    }
}


//int flag = 0;
/* sample a x in [m] with weight w w indices: 0 ~ (m-1). */
int SampleXWeight(double weight[], int m){
    
    int x;
    double sum=0.0;
    for (int i = 0; i < m; i++) {
        sum += weight[i];
    }
    // cout<<"After normalization "<<endl;
    for (int i = 0; i < m; i++) {
        weight[i] = weight[i]/sum;
        //   cout<<weight[i]<<' '<<endl;
    }
    //  cout<<endl;
    
    x = DirectSimulate(weight);
    //  cout<<x<<endl;
    return x;
}


/*factorial of n = n * (n-1) * ... * 2 * 1 */
int Factorial(int n){
    int ans=1;
    for (int i = 2; i <= n; i++) {
        ans *= i;
    }
    return ans;
}

/*qsort for x*/
int compare(const void *a, const void *b){
    return (*(int*)a - *(int*)b);
}
/*sum_{sigma in S_A} prod_{1\leq j \leq ta} [x_{sigma_i = x_B_i}]=Per (SM)*/
//void ConstructBipar(int x[], int a[], int b[], int t, int SM[][MAXY]){
//    
//    for (int i = 0; i < t; i ++) {
//        for (int j = 0; j < t; j ++) {
//            if (x[a[i]] == x[b[j]]) {
//                SM[i][j] = 1;
//            }
//            else SM[i][j] = 0;
//        }
//    }
//    
//}

/* Compute the number of perfect matchings between xa and xb.*/
double BiparPer(int x[], int a[], int b[], int ta){
	if(ta == 0)return 1;
	int ans;
	int y[100], z[100];
	
	for(int i = 0; i < ta; i ++){
		y[i] = x[a[i]];
		z[i] = x[b[i]];
	} 
	qsort(y,ta,sizeof(int),compare);
	qsort(z, ta, sizeof(int), compare);
	ans = 1;
	int temp, count[100];
	temp = y[0];
	int k = 0;
	count[0] = 0;
	for(int i = 0; i < ta; i++){
		if(y[i] != z[i]){ /*There exists a element u in xa and not in xb, then per = 0 */
			return 0;
		}
		if(y[i] == temp){  /* counting the outcoming times of i in xa. Per = the factorial of which. */
			count[k]++;
		}
		else {
			temp = y[i];
			k++;
			count[k]=1;
		}
	}
	for(int i = 0; i < k + 1; i ++){
		ans *= Factorial(count[i]);
	}
	return ans;
}


//map<string, Complex > mapHaf; // mapHaf<index, hafnian>


/*check of the correctness.*/
void GBS(int n, int m, int x[]){
    
    double weight[MAXX];
    Complex z;
    Complex hafz;
    Complex za,zb,zea,zeb;
    
  //  FILE *fp;
  //  fp = fopen("Weight.txt", "w+");
    // Complex ztemp;
    for (int i = 0; i < m; i++) {
        weight[i]=1.0;
    }
    int binomial[50][50];
    for (int j = 0 ; j < n/2; j ++) {
        binomial[j][0] = 1;
        for (int k = 1; k <= n/2; k ++) {
            binomial[j][k] = 0;
        }
    }
    for (int j = 0 ; j <= n/2; j ++) {
        binomial[j][0] = 1; /* binomial[ j, k ] = C(m/2 - 1 + j, k )  */
        for (int k = 1; k <= j ; k ++) {
            binomial[j][k] = (binomial[j][k - 1] * (m/2 + j - k ) )/ k;
        }
    }
    // compute C(m/2 - 1 + j, k) for 0 <= j <=n/2 and 0 <= k <= j
    
    x[0]=SampleXWeight(weight,m);
   // x[0] = 0;
    
    for (int k = 1; k < n; k ++) {
        
        for (int l = 0; l < m; l ++){
            
            Complex zsx;
            
            x[k] = l;
            weight[l] = 0;
            
            hafz.real = 0;
            hafz.image = 0;
            
            //   #pragma omp parallel for
            for (int ja = max(0, 2 * (k + 1) - n); ja <= k + 1; ja += 2) {
                for (int jb = max(0, 2 * (k + 1) - n); jb <= ja; jb += 2) {
                    
                    string bitsa(ja, 1);
                    bitsa.resize(k + 1, 0);
                    
                    int sa[MAXX],sb[MAXX], sua[MAXX],sub[MAXX],sea[MAXX],seb[MAXX], sta[MAXX], stb[MAXX], SW[MAXY][MAXY];
                    int counta, countb, countea, counteb, countta, counttb,counttea, countteb, temp;
                    
                    // a = max((k + 1) % 2, 3 * (k + 1) - (ja + jb) - n) ;
                    za.real = za.image = 0;
                    
                    //   #pragma omp parallel do
                    do{
                        string sx;
            //            map< string, Complex>::iterator it;
                        counta = 0;
                        countta = 0;
                        
                        for (int i = 0; i < k + 1; i ++) {
                            if (bitsa[i]) {
                                sa[counta ++] = i;
                                
                            }
                            else{
                                sua[countta ++] = i;
                            }
                            
                        }
                        
                        string bitsb(jb, 1);
                        bitsb.resize(k + 1, 0);
                        
                        zb.real = zb.image = 0;
                        
                        // #pragma omp parallel do
                        do{
                            countb = 0;
                            counttb = 0;
                            for (int i = 0; i < k + 1; i ++) {
                                if (bitsb[i]) {
                                    sb[countb ++] = i;
                                    
                                }
                                else {
                                    sub[counttb ++] = i;
                                }
                                
                            }
                            Complex sumz;
                            sumz.real = 0;
                            sumz.image = 0;
                            for (int ta = max( (k + 1) % 2, 3 * (k + 1) - (ja + jb) - n); ta <= k + 1 - max(ja,jb); ta += 2) {
                                
                                string bitsfa(ta, 1);
                                bitsfa.resize(k + 1 - ja, 0);
                                
                                //int factor2 = Binom((n - (k + 1) + ta + m)/2 - 1, k  + (m - ja - jb)/2 );
                                int factor;
                                int upper = (n - (k + 1) + ta + m)/2 - 1;
                                int lower = min(k + (m - ja - jb)/2 , upper - k  - (m - ja - jb)/2);
                                factor = binomial[upper - m/2 + 1][lower];
                                //                                if(factor != factor2){
                                //                                    cout<<"Factor error!!"<<upper<<" and "<<lower <<endl;
                                //
                                //                                }
                                zea.real = zea.image = 0;
                                do{
                                    countea = 0;
                                    counttea = 0;
                                    
                                    for (int i = 0; i < k + 1 - ja; i ++) {
                                        if (bitsfa[i]) {
                                            sta[countea ++] = sua[i];
                                        }
                                        else{
                                            sea[counttea ++] = sua[i];
                                        }
                                    }
                                    string bitsfb(ta, 1);
                                    bitsfb.resize(k + 1 - jb, 0);
                                    zeb.real = zeb.image = 0;
                                    do{
                                        counteb = 0;
                                        countteb = 0;
                                        for (int i = 0; i < k + 1 - jb; i ++) {
                                            if(bitsfb[i]){
                                                stb[counteb ++] = sub[i];
                                            }
                                            else {
                                                seb[countteb ++] = sub[i];
                                            }
                                        }
                                       // ConstructBipar(x, sta, stb, ta, SW);
									//	temp = Permanent(SW, ta);
										temp = BiparPer(x,sta, stb, ta);
                                        if(temp == 0){
                                        	z.real = 0;
                                        	z.image = 0;
										}
										else z = Hafnian(W, x, seb, k + 1 - jb - ta);
                                        
                                        z.real = z.real * temp;
                                        z.image = z.image * temp;
                                        zeb.real += z.real;
                                        zeb.image -= z.image;
                                    }while (prev_permutation(bitsfb.begin(), bitsfb.end()));
                                    
                                    z = Hafnian(W, x, sea, k + 1 - ja - ta);
                                    z = MultiComplex(z, zeb);
                                    zea.real += z.real;
                                    zea.image += z.image;
                                    
                                }while (prev_permutation(bitsfa.begin(), bitsfa.end()));
                                zea.real *= factor;
                                zea.image *= factor;
                                
                                sumz.real += zea.real;
                                sumz.image += zea.image;
                                
                            }//end-for
                            z = Hafnian(W, x, sb, jb);
                            
                            z.image  = -1 * z.image;
                            z = MultiComplex(z, sumz);
                            
                            zb.real += z.real;
                            zb.image += z.image;
                            
                        }while (prev_permutation(bitsb.begin(), bitsb.end()));
                        z = Hafnian(W, x, sa, ja);
                        z = MultiComplex(z, zb);
                        
                        za.real += z.real;
                        za.image += z.image;
                        
                    }while (prev_permutation(bitsa.begin(),bitsa.end()));
                    
                    hafz.real += za.real;
                    hafz.image += za.image;
                    
                    if (jb < ja) {
                        hafz.real *= 2;
                    }
                }//end-jb
            }//end-ja
            weight[l] = hafz.real;
        //    fprintf(fp, "%lf ", weight[l]);
        }//end -for-l
        x[k] = SampleXWeight(weight, m);
      //  x[k] = k;
      //  fprintf(fp, "\n");
    }//end-for-k
    
   // fclose(fp);
}


int main() {
    // insert code here...
    int n,m;
    int x[MAXX];
    //double A[MAXX][MAXX];
    
    double **A = new double *[MAXW];
    for (int j = 0; j < MAXW; j ++){
        A[j] = new double [MAXW * 2];
    }
    W = new double *[MAXW];
    for (int j = 0; j < MAXW; j ++){
        W[j] = new double[MAXW * 2];
    }
    printf("Please input the size of unitary matrix m: (We will get a m*m matrix.) \n");
    cin>>m;
    printf("Please input the size of photons \n");
    cin>>n;
   	SchmidtOrth(m, A);
	SymmetricM(A, m);
 //   FILE * fp;
//    fp = fopen("timeSec.txt", "w+");

    clock_t t = clock();
    GBS(n,m,x);
    qsort(x, n, sizeof(int), compare);
    t = clock() - t;
    cout<<"The sample obtained by classical simulation of GBS."<<endl;
    for (int i = 0; i < n; i++) {
        cout<<x[i]<<" ";
    }
    cout<<endl;
    
    printf("Time : %f\n",((float)t)/CLOCKS_PER_SEC);
//    int T = 300;
//    
//    clock_t start, end;
//   
//    for(int n = 2; n <= 12; n += 2){
//
//	m = n * n;
////	cout<<m<<endl;
//	SchmidtOrth(m, A);
//	SymmetricM(A, m);
//	start = time(NULL);
//	t = clock();
//        for(int i = 0; i < T; i ++){
//		GBS(n, m, x);
//	}
//	t = clock() - t;
//	end = time(NULL);
////	cout<<n<<": "<<difftime(end, start)<<endl;
//	cout<<n<<", "<<m<<": "<<((float)t)/CLOCKS_PER_SEC/T<<endl;
//	fprintf(fp, "%f\n", ((float)t)/CLOCKS_PER_SEC/T);
//    }
    
   // cout<<"The number of photons are "<<n<<endl;
   // for (int i = 0; i < n; i++) {
   //     cout<<x[i]<<" ";
   // }
   // cout<<endl;
//   fclose(fp);
    for (int j = 0; j < MAXW; j ++){
        delete [] W[j];
    }
    for (int j = 0; j <MAXW; j ++){
        delete [] A[j];
    }
    return 0;
}
