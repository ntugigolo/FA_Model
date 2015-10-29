#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <time.h>
#include <random>
#include <algorithm>
#include <chrono>
#define pi 3.14159265359
using namespace std;

double max(double a1,double b1){
	if(a1>b1) return a1;
	else return b1; 
}

double uniformRandom(){
  return (double)(rand())/(double)(RAND_MAX);
}

void mo_callprice_eur(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r,const double& e,const double& m){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	double rand ; 
	double call = 0;
	double call_1 = 0;
	double call_2 = 0;
	double price,multipri,price1,price2;
	double var1,var2,x;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rng(seed);
	//std::random_device rd;
	//std::mt19937 rng;
	std::normal_distribution<> ND(0.,1.);
	for(int i = 1;i<= m;i++){
		price = S;//multipri = S;
		price1 = S + e;
		price2 = S - e;
		for(int j = 1;j <= n;j++){
		double phi = ND(rng);
		price = price*exp((r-pow(s,2)/2)*(t/n)+s*sqrt(t/n)*phi);
		price1 = price1*exp((r-pow(s,2)/2)*(t/n)+s*sqrt(t/n)*phi);
		price2 = price2*exp((r-pow(s,2)/2)*(t/n)+s*sqrt(t/n)*phi);
		}
		call = price >= X ? call+1 : call;
		/*double com = price-X;
		call = call + std::max(com,0.0);*/
		call_1 = price1 >= X ? call_1+1 : call_1;
		call_2 = price2 >= X ? call_2+1 : call_2;
		} 
	cout<<"The European call price in Monte Carlo:"<< call/exp(r*t)/m <<endl;
	cout<<"Delta in Monte Carlo:"<< ((call_1/exp(r*t)/m)-(call_2/exp(r*t)/m))/(2*e) <<endl;	
	cout<<"Gamma in Monte Carlo:"<< ((call_1/exp(r*t)/m)-2*(call/exp(r*t)/m) +(call_2/exp(r*t)/m))/(pow(e,2)) <<endl;	
}

void tri_callprice_eur(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r,const double& e){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	int h = log(S/X)/(s*sqrt(t/n)); 
	//cout<<"sDG"<<h<<endl;
	// h = 5;
	double landa = sqrt(pi/2);
	//double landa = sqrt(2);
//	double landa = log(S/X)/(h*s*sqrt(t/n));
	 //landa = 1.17797;
	u = exp(landa*s*sqrt(t/n));
	//u = 1.0213;
	d= 1/u;
	double uu=u*u;
	double p_up = 1/(2*pow(landa,2))+(r+pow(s,2))*sqrt(t/n)/(2*landa*s);
	double p_down= 1/(2*pow(landa,2))-(r-2*pow(s,2))*sqrt(t/n)/(2*landa*s);
	double p_m = 1-p_up-p_down;
	//cout<<"D:"<<p_m<<endl;
	vector<double> prices(n*2+1);
	prices[0]=S*pow(u,n);
	for(int i=1;i<=n*2;++i) prices[i] = d*prices[i-1];
	//cout<<"dsf:"<<prices[2000]<<endl;
	vector<double> call_values(n*2+1);
	for(int i=0;i<=n*2;++i){
	 if(prices[i]>= X) call_values[i] = 1;
	 else call_values[i] = 0; 		
	 }
	 vector<double> call_gamma;
	for(int step=n-1; step>=0; --step){
		for(int i=0;i<=step*2;++i){
			call_values[i] = (p_up*call_values[i]+p_down*call_values[i+2]+p_m*call_values[i+1])*Rinv;
			}
			if(step == 1){
				for(int i = 0;i<=4;i++){
				call_gamma.push_back(call_values[i]);	
				}
			}
	}
	cout<<"The European call price in trinomial tree:"<< call_values[0]<<endl;
	cout<<"Delta in trinomial tree:"<< (call_gamma[0]-call_gamma[2])/(S*u-S*d)<<endl;
	cout<<"Gamma in trinomial tree:"<< ((call_gamma[0]-call_gamma[1])/(S*u-S)-(call_gamma[1]-call_gamma[2])/(S-S*d))/((S*u-S*d)/2)<<endl;
}
int main(){
	double S,X,s,t,n,r,m,e;
	cout<<"stock price at time = 0 :";
	cin>>S;
	cout<<"strike price :";
	cin>>X;
	cout<<"annual volatility in percentage :";
	cin>>s;
	cout<<"maturity in years :";
	cin>>t;
	cout<<"the number of periods :";
	cin>>n;
	cout<<"interest rate in percentage :";
	cin>>r;
	cout<<"number of sample paths of simulation(for Monte Carlo):";
	cin>>m;
	cout<<"e (for Monte Carlo):";
	cin>>e;
	/*S=100;
	X=90;
	s=40;
	r=5;
	t=2;
	n=1000;
	m = 1000000;
	e = 0.001;*/
	s=s/100;
	r=r/100;
	tri_callprice_eur(S,X,s,t,n,r,e);
	//mo_callprice_eur(S,X,s,t,n,r,e,m);
	//system("pause");
	return 0;	
}

