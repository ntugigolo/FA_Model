#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <cstdio>
#include <math.h>

using namespace std;

class FA4{
	public:
		FA4();
		~FA4();
		double max(double,double);
		double prob(double,double,double);
		void callprice_eur(const double&,const double&,const double&,const double&,const int&,const double&,const int&,const double&,const int&);
};

FA4::FA4(){}

FA4::~FA4(){}

double FA4::max(double a1,double b1){
	if(a1>b1) return a1;
	else return b1; 
}
double FA4::prob(double R,double u,double d){
	return (R-d)/(u-d);
}

void FA4::callprice_eur(const double& S,const double& X,const double& s,const double& t,const int& n,const double& r,const int& k,const double& H,const int& Delta){
	double R,u,d,C;
	R = exp(r*(t/n));
	double Rinv=1.0/R;
	u = exp(s*sqrt(t/n));
	d= 1.0/u;
	double prob_up = prob(R,u,d);
	double prob_down=1-prob_up;	
	double eps = 0.0001;

	vector<double> price_states_end(k*(n+1));
	vector<double> price_max_end(n+1);
	vector<double> price_min_end(n+1);
	for(int j=n;j>=0;j--){
		price_max_end[n-j]=(S*(1.0-pow(u,n-j+1.0))/(1.0-u) + S*pow(u,n-j)*d*(1.0-pow(d,j))/(1.0-d))/(n+1.0);	
		price_min_end[n-j]=(S*(1.0-pow(d,j+1.0))/(1.0-d) + S*pow(d,j)*u*(1.0-pow(u,n-j))/(1.0-u))/(n+1.0);
		for(int i=0;i<k;i++){	
		price_states_end[k*(n-j)+i] = i/(k-1.0)*price_max_end[n-j] + (k-i-1.0)/(k-1.0)*price_min_end[n-j];	
		}
	}

	vector<double> call_values_end(k*(n+1.0));
	for(int i=0;i< k*(n+1.0);++i) {
		if(price_states_end[i] >= H){
			call_values_end[i] = 0;
		}
		else{
		call_values_end[i] = max(0,(price_states_end[i]-X));
	}
	}
	double p_up,p_down,c_up,c_down;
	vector<double> price_states_pre;
	vector<double> call_values_pre;	
	vector<double> call_values_p1;

	for(int step=n-1; step>=0; --step){
		vector<double> price_states(k*(step+1.0));
		vector<double> price_max(step+1.0);
		vector<double> price_min(step+1.0);
		vector<double> call_values(k*(step+1.0));
		
		for(int j=step;j>=0;j--){
		price_max[step-j]=(S*(1.0-pow(u,step-j+1.0))/(1.0-u) + S*pow(u,step-j)*d*(1.0-pow(d,j))/(1.0-d))/(step+1.0);	
		price_min[step-j]=(S*(1.0-pow(d,j+1.0))/(1-d) + S*pow(d,j)*u*(1.0-pow(u,step-j))/(1.0-u))/(step+1.0);
		
		for(int i=0;i<k;i++){
		price_states[k*(step-j)+i] = i/(k-1.0)*price_max[step-j] + (k-i-1.0)/(k-1.0)*price_min[step-j];	

		if(step==n-1.0){
		p_up = ((step+1.0)*price_states[k*(step-j)+i]+S*pow(u,step+1-j)*pow(d,j))/(step+2.0);
		p_down = ((step+1.0)*price_states[k*(step-j)+i]+S*pow(u,step-j)*pow(d,j+1))/(step+2.0);
			for(int m=0;m<k;m++){
				if( fabs(p_down - price_states_end[k*(step-j)+m] ) < eps ){
					c_down = call_values_end[k*(step-j)+m];	
					break;
				}
				if( fabs(p_down - price_states_end[k*(step-j)+m] ) > eps && p_down < price_states_end[k*(step-j)+m] ){
					double w_d = ((price_states_end[k*(step-j)+m]-p_down)/(price_states_end[k*(step-j)+m]-price_states_end[k*(step-j)+m-1.0]));
					if( p_down >= H) c_down = 0;
					else	c_down = w_d*call_values_end[k*(step-j)+m-1.0] + (1.0-w_d)*call_values_end[k*(step-j)+m];
						break;
						}	  
						}
			for(int m=0;m<k;m++){			
				if( fabs(p_up - price_states_end[k*(step-j+1)+m] )< eps){
					c_up = call_values_end[k*(step-j+1.0)+m];
					break;	
					}
				if( fabs(p_up - price_states_end[k*(step-j+1)+m] ) > eps && p_up < price_states_end[k*(step-j+1)+m]){
					double w_u = ((price_states_end[k*(step-j+1.0)+m]-p_up)/(price_states_end[k*(step-j+1.0)+m]-price_states_end[k*(step-j+1.0)+m-1.0]));
					if( p_up >= H) c_up = 0;
					else c_up = w_u*call_values_end[k*(step-j+1.0)+m-1.0] + (1.0-w_u)*call_values_end[k*(step-j+1.0)+m];
				break;
						}
					}
		if(price_states[k*(step-j)+i] >=  H)  call_values[k*(step-j)+i] = 0;		
		call_values[k*(step-j)+i] =  (prob_up*c_up + prob_down*c_down)/R;
			}	
		else{
		p_up = ((step+1.0)*price_states[k*(step-j)+i]+S*pow(u,step+1.0-j)*pow(d,j))/(step+2.0);
		p_down = ((step+1.0)*price_states[k*(step-j)+i]+S*pow(u,step-j)*pow(d,j+1.0))/(step+2.0);
		
			for(int m=0;m<k;m++){
				if( fabs(p_down - price_states_pre[k*(step-j)+m] )< eps){	
					c_down = call_values_pre[k*(step-j)+m];	
					break;
				}
				if( fabs(p_down - price_states_pre[k*(step-j)+m] ) > eps && p_down < price_states_pre[k*(step-j)+m]){
					double w_d = ((price_states_pre[k*(step-j)+m]-p_down)/(price_states_pre[k*(step-j)+m]-price_states_pre[k*(step-j)+m-1.0]));
					if (p_down >= H) c_down = 0;
					else c_down = w_d*call_values_pre[k*(step-j)+m-1.0] + (1.0-w_d)*call_values_pre[k*(step-j)+m];
					break;
						}
						}
			for(int m=0;m<k;m++){			
			if( fabs(p_up - price_states_pre[k*(step-j+1)+m] )< eps){
					c_up = call_values_pre[k*(step-j+1.0)+m];	
					break;
					}
			if( fabs(p_up - price_states_pre[k*(step-j+1)+m] ) > eps && p_up < price_states_pre[k*(step-j+1)+m]){
					double w_u = ((price_states_pre[k*(step-j+1.0)+m]-p_up)/(price_states_pre[k*(step-j+1.0)+m]-price_states_pre[k*(step-j+1.0)+m-1.0]));
					if (p_up >= H) c_up = 0;
					else c_up = w_u*call_values_pre[k*(step-j+1.0)+m-1.0] + (1.0-w_u)*call_values_pre[k*(step-j+1.0)+m];	
					break;
						}
				}
		if(price_states[k*(step-j)+i] >=  H) call_values[k*(step-j)+i] = 0;		
		call_values[k*(step-j)+i] =  (prob_up*c_up + prob_down*c_down)/R;	
			}			
			}	
		}
		if(step == 1){
		call_values_p1 = call_values;	
		}
		price_states_pre = price_states;
		call_values_pre = call_values;
		}
	cout<<"The European call price :"<<call_values_pre[0]<<endl;
	cout<<"Delta :"<< (call_values_p1[k]-call_values_p1[0])/(S*u-S*d);
	} 
int main(){
	double S,X,s,t,r,H,delta;
	int k,n;
	FA4 A;
	cout<<"stock price at time = 0 :";
	cin>>S;
	cout<<"strike price :";
	cin>>X;
	cout<<"Barrier :";
	cin>>H;
	cout<<"maturity in years :";
	cin>>t;
	cout<<"annual volatility in percentage :";
	cin>>s;
	cout<<"continuously compounded annual interest rate :";
	cin>>r;
	cout<<"number of periods :";
	cin>>n;
	cout<<"number of states per node :";
	cin>>k;

	s=s/100;
	r=r/100;
	delta =1;
	A.callprice_eur(S,X,s,t,n,r,k,H,delta);
	system("pause");
	return 0;
}

