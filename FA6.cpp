#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>
using namespace std;
class FA6{
	public:
		FA6();
		FA6(double&,double&,double&,double&,double&,double&,double&,double&,double&,double&);
		~FA6();
		double max(double,double);
		double spot_rate(const double&);
		double delta();
		double con_price(const double&,const double& ,const double&);
		double short_rate(const double&);
		//double callprice_eur(const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&);
		double callprice_eur_con(const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&);
	private:
		double X,s,x,n,y,E,F,G,prob,del_t;
};

FA6::~FA6(){}

FA6::FA6(double& x_,double& y_,double& E_,double& F_,double& G_,double& s_,double& n_,double& X_,double& select,double& prob_){
	del_t = y_/n_;
	x=x_; y=y_; E = E_; F = F_ ; G = G_ ; s = s_ ; n = n_ ; X = X_ ; prob = prob_;
	//if (select == 1) callprice_eur(x,y,E,F,G,s,n,X);
	if (select == 2) callprice_eur_con(x,y,E,F,G,s,n,X);
}

double FA6::max(double a1,double b1){
	if(a1>b1) return a1;
	else return b1; 
}
double FA6::spot_rate(const double& t){
	double spot = E - F*exp(-1*G*t);
	return spot;
}
double FA6::delta(){
	//cout<<"s:"<<s<<" del_t:"<<del_t<<"prob:"<<prob<<endl;
	double del = exp((-1*s*pow(del_t,1.0/2.0))/sqrt(prob*(1-prob)));
	return del;
}
/*double FA6::callprice_eur(const double& x,const double& y,const double& E,const double& F,const double& G,const double& s,const double& n,const double& X){
	vector<double> prices(n+1);
	//vector<double> Rate(n+1);
	vector<double> sub_Rate(n+1);
	vector<double> call(n+1);
	prices[0] = 1/pow((1+spot_rate(E,F,G,y)),n);
	double P_T = prices[0];
	double P_t = 1/pow((1+spot_rate(E,F,G,x)),n/y*x);
	double prob = 0.5;
	for (int i =0;i<=n/y*x;i++){
		prices[i] = P_T/P_t	* pow(delta(s,y/n,prob),(y-x)*(n/y*x-i));
		for (int j=0;j<n/y*x;j++){
		prices[i] *= (prob + (1-prob) * pow(delta(s,y/n,prob),x-(j+1)*y/n) )/(prob + (1-prob) * pow(delta(s,y/n,prob),y-(j+1)*y/n));
		}
	}
	for (int i =0;i<=n/y*x;i++){
	call[i] = max(prices[i]-X,0);
	}
	double P_dt,P_ut;
	vector<double> sub_prices(n+1);
	for(int k = n/y*x-1;k>=0;k--){
		 P_dt = 1/pow((1+spot_rate(E,F,G,k*y/n)),k);
		 P_ut = 1/pow((1+spot_rate(E,F,G,(k+1)*y/n)),k+1);
		for (int i =0;i<=k;i++){
		sub_prices[i] = P_ut/P_dt*pow(delta(s,y/n,prob),(y/n)*(k-i));
		for (int j=0;j<=k;j++){
		sub_prices[i] *= (prob + (1-prob) * pow(delta(s,y/n,prob),k*y/n-(j+1)*y/n) )/(prob + (1-prob) * pow(delta(s,y/n,prob),(k+1)*y/n-(j+1)*y/n));
		}
		sub_Rate[i] = -1*log(sub_prices[i])*n/y;
		call[i] = (call[i+1]*prob+call[i]*(1-prob))/(1+sub_Rate[i]);
		} 		
	}
	cout<<"The European call option on a zero-coupon bond :";
		cout<< call[0]<<endl;
}*/
double FA6::con_price(const double& state,const double& t,const double& T){
	double P_T = 1/exp(T*spot_rate(T));
	double P_t = 1/exp(t*spot_rate(t)); 
	/*double P_T = 1/exp(T/del_t*spot_rate(T));
	double P_t = 1/exp(t/del_t*spot_rate(t));*/
	double prices;
		prices = P_T/P_t*pow(delta(),(T-t)*(t/del_t-state));
		for (int j=0;j<t/del_t;j++){
		prices *= (prob + (1-prob) * pow(delta(),t-(j+1)*del_t) )/(prob + (1-prob) * pow(delta(),T-(j+1)*del_t));
		}
	return prices;	
}
double FA6::callprice_eur_con(const double& x,const double& y,const double& E,const double& F,const double& G,const double& s,const double& n,const double& X){
	vector<double> prices(n+1);
	vector<double> sub_Rate(n+1);
	vector<double> call(n+1);
	for (int i = 0;i <= n/y*x;i++){
		prices[i] = con_price(i,x,y);
	}
	for (int i =0;i<=n/y*x;i++){
	call[i] = max(prices[i]-X,0);
	}
	vector<double> sub_prices(n+1);
	for(int k = x/del_t-1;k>=0;k--){
		for (int i =0;i<=k;i++){
		sub_prices[i] = con_price(i,k*del_t,(k+1)*del_t);
		sub_Rate[i] = short_rate(sub_prices[i]);
		call[i] = (call[i+1]*prob+call[i]*(1-prob))/exp(del_t*sub_Rate[i]);
		} 		
	}
		cout<<"The European call option on a zero-coupon bond :";
		cout<<call[0]*100<<" (% of par)";
}

double FA6::short_rate(const double& sub_pri){
	double rate = -1*log(sub_pri)/del_t;
	return rate;
}

int main(){
	double X,s,x,n,y,E,F,G,prob;
	cout<<"year(call option) :";
	cin>>x;
	cout<<"year(bond) :";
	cin>>y;
	cout<<"E :";
	cin>>E;
	cout<<"F :";
	cin>>F;
	cout<<"G :";
	cin>>G;
	cout<<"constant annualized volatility of short rate :";
	cin>>s;
	cout<<"The number of time periods of the tree :";
	cin>>n;
	cout<<"Strike price in % of par :";
	cin>>X;
	prob = 0.5;
	double mode = 2;//1 for discrete, 2 for continous;
/*	x = 1;
	y = 2, E = 0.08, F = 0.05, G = 0.18, s = 10, n = 30,X = 90;*/
	s=s/100;
	X=X/100;
	FA6 Ans(x,y,E,F,G,s,n,X,mode,prob);
	//callprice_eur_con(x,y,E,F,G,s,n,X,mode,prob);
	system("pause");
	return 0;
}

