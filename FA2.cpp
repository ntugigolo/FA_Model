#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>

using namespace std;
int factorial(double j){
	if(j==0)
	return 1;
	else
	return j*factorial(j-1);	
}
double max(double a1,double b1){
	if(a1>b1) return a1;
	else return b1; 
}
double prob(double R,double u,double d){
	return (R-d)/(u-d);
}
double callprice_eur(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;
	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n);
	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> call_values(n+1);
	for(int i=0;i<=n;++i) call_values[i] = max(0,(prices[i]-X));
	for(int step=n-1; step>=0; --step){
		for(int i=0;i<=step;++i){
			call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
		}
	}/*
	for(int j=0;j<=n;j++){
	C += 1/pow(R,n)*factorial(n)/factorial(j)/factorial(n-j)*pow(p,j)*pow(1-p,n-j)*max(0,pow(u,j)*pow(d,n-j)*S-X);		
	}*/
	//return C;
	return call_values[0];
}
double putprice_eur(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;
	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n);
	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> put_values(n+1);
	for(int i=0;i<=n;++i) put_values[i] = max(0,(X-prices[i]));
	for(int step=n-1; step>=0; --step){
		for(int i=0;i<=step;++i){
			put_values[i] = (p_up*put_values[i+1]+p_down*put_values[i])*Rinv;
		}
	}
	return put_values[0];
}
double callprice_ama(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;
	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n);
	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> call_values(n+1);
	for(int i=0;i<=n;++i) call_values[i] = max(0,(prices[i]-X));
	for(int step=n-1; step>=0; --step){
		for(int i=0;i<=step;++i){
			call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
			prices[i] = d*prices[i+1];
			call_values[i] = max(call_values[i],prices[i]-X);
		}
	}
		return call_values[0];
}
double putprice_ama(const double& S,const double& X,const double& s,const double& t,const double& n,const double& r){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1/R;
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;

	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n);
	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> put_values(n+1);
	for(int i=0;i<=n;++i) put_values[i] = max(0,(X-prices[i]));
	for(int step=n-1; step>=0; --step){
		for(int i=0;i<=step;++i){
			put_values[i] = (p_up*put_values[i+1]+p_down*put_values[i])*Rinv;
			prices[i] = d*prices[i+1];
			put_values[i] = max(put_values[i],X-prices[i]);
		}
	}
		return put_values[0];
}

int main(){
	double S,X,s,t,n,r;
	
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
	s=s/100;
	r=r/100;
	cout<<"The European call price :"<<callprice_eur(S,X,s,t,n,r)<<endl;
	cout<<"The European put price :"<<putprice_eur(S,X,s,t,n,r)<<endl;
	cout<<"The Amarican call price :"<<callprice_ama(S,X,s,t,n,r)<<endl;
	cout<<"The Amarican put price :"<<putprice_ama(S,X,s,t,n,r)<<endl;
	system("pause");
	return 0;
	
}

