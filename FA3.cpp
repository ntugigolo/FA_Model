#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include<algorithm>
using namespace std;

int WeekDate(int year,int month,int day){
	int w;
	int m;
	if(month == (1||2)){ m = month +12; year-=1;}
	else m=month;
	
	int c = year/100;
	int y=year%100;
	
	w = (y + floor(y/4) + floor(c/4) - 2*c + floor(26*(m+1)/10) +day-1);
	w = ( w % 7 + 7 ) % 7;
	return w;
}

void StringToDate(string date,int& year,int& month,int& day) { 
	year = atoi((date.substr(0,4)).c_str()); 
	month = atoi((date.substr(5,2)).c_str()); 
	day = atoi((date.substr(8,2)).c_str());}

bool IsLeap(int year) { 
 return
  (year%4 == 0) && (year%100 != 0) || (year%400 == 0); }

int DayInYear(int year,int month,int day) { 

	int DAY[12] = {31,28,31,30,31,30,31,31,30,31,30,31}; 
	if (IsLeap(year)) 
		DAY[1] = 29; 
	for(int i=0;i<month-1;++i) 
		day += DAY[i]; 
	return day; }

int DaysBetween2Date(string date1,string date2){ 
	int year1,month1,day1; 
	int year2,month2,day2;
	StringToDate(date1,year1,month1,day1);
	StringToDate(date2,year2,month2,day2);

	if (year1==year2 && month1==month2) { return day1>day2 ? day1-day2 : day2-day1; }
	else if (year1 == year2) { 
	int d1,d2; 
	d1 = DayInYear(year1,month1,day1);
	d2 = DayInYear(year2,month2,day2); 
	return d1>d2 ? d1-d2 : d2-d1; }
else { 
if(year1 > year2){ 
swap(year1,year2); 
swap(month1,month2); 
swap(day1,day2); }

int d1,d2,d3; 
if(IsLeap(year1)) 
	d1 = 366 - DayInYear(year1,month1,day1); 
else 
	d1 = 365 - DayInYear(year1,month1,day1); 
d2 = DayInYear(year2,month2,day2); 
d3 = 0;

for(int year = year1+1; year < year2; year++) {
 
if (IsLeap(year)) 
	d3 += 366; 
else 
	d3 += 365; 
} 
return d1 + d2 +d3; 
}
}
double prob(double R,double u,double d){
	return (R-d)/(u-d);
}
double callprice_Bermu(const double& S,const double& X,const double& s,string& T0,string& T1,const int& m,const double& r,string& T2,string& T3){
	
	double R,u,d,C;
	int year0,month0,day0; 
	int year2,month2,day2;
	int year1,month1,day1; 
	int year3,month3,day3;
	StringToDate(T0,year0,month0,day0);
	StringToDate(T1,year1,month1,day1);
	StringToDate(T2,year2,month2,day2);
	StringToDate(T3,year3,month3,day3);
	int w0 = WeekDate(year0,month0,day0);
	int w1 = WeekDate(year1,month1,day1);
	int w2 = WeekDate(year2,month2,day2);
	int w3 = WeekDate(year3,month3,day3);
	int Total_days = DaysBetween2Date(T0,T1);
	int Total_days_T2 = DaysBetween2Date(T0,T2);
	int Total_days_T3 = DaysBetween2Date(T0,T3);
	
	vector<int> Weekend;
	int F = (6-w0)*m;
	int Lim;
	if(w0 < w1 ){
		 Lim =  (Total_days/7);
	}
	else{
		 Lim = (Total_days/7 +1 );
	}
	
	for(int i=0;i< Lim;i++){
		
		Weekend.push_back(F+1);
		Weekend.push_back(F+2);
		Weekend.push_back(F+3);
		Weekend.push_back(F+4);
		F += 7*2;
	}
	/*for(vector<int>::iterator it=Weekend.begin();it!=Weekend.end();it++){
		
		cout<<"Sephiria The weekend:"<<*it<<endl;
	}*/
	int Trade_days;	
	if(w0 < w1 ){
		Trade_days =Total_days - (Total_days/7)*2 + 1;
	}
	else{
		Trade_days =Total_days - (Total_days/7 +1 )*2 + 1;
	}
	double First,End;
	if(w0 < w2 ){
		First =Total_days_T2 - (Total_days_T2/7)*2;
	}
	else{
		First =Total_days_T2 - (Total_days_T2/7 + 1 )*2;
	}
	if(w0 < w3 ){
		End =Total_days_T3 - (Total_days_T3/7)*2 + 1;
	}
	else{
		End =Total_days_T3 - (Total_days_T3/7 + 1 )*2 + 1;
	}
	int holiday;
	if(w0 < w2 ){
		holiday =(Total_days/7)*2;
	}
	else{
		holiday =(Total_days/7 + 1 )*2;
	}	
	Total_days += 1;
	Total_days_T3 += 1;
	double t2 =Total_days/365.0;
	double n2 = Total_days*m;
	/*double First,End;
	First =Total_days_T2;
	 End =Total_days_T3;*/
	//cout<<"Total_days:"<<End<<endl;
	double t =Trade_days/261.0;
	int n = Trade_days*m;
	/*double t =44.0/261.0;
	double n = 44.0*m;*/
	//	cout<<"GDFG:"<<n<<endl;
	First *=m;
	End*=m;
//cout<< "t"<<setprecision(6)<<t<<endl;
	R = exp(r*(t2/n2));
	double Rinv=1/R;
	
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;
	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n)*pow(R,holiday*m);
	
	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> call_values(n+1);
	int check = 0;
	int k =0;
	int ite=0;
	for(int i=0;i<=n;++i) call_values[i] = max(0.0,(prices[i]-X));
	for(int step=n-1; step>=0; --step){

		for(int i=0;i<=step;++i){
			call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
			prices[i] = d*prices[i+1];
		if(step > First&&step <= End)	call_values[i] = max(call_values[i],prices[i]-X);
			if(((n-w1*m)-step) % 10 == 0){
				call_values[i] = call_values[i]*pow(Rinv,2*m);
				prices[i] = (prices[i])*pow(Rinv,2*m);
		}
	}
	}
		return call_values[0];
}

double putprice_Bermu(const double& S,const double& X,const double& s,string& T0,string& T1,const int& m,const double& r,string& T2,string& T3){
	double R,u,d,C;
	int year0,month0,day0; 
	int year2,month2,day2;
	int year1,month1,day1; 
	int year3,month3,day3;
	StringToDate(T0,year0,month0,day0);
	StringToDate(T1,year1,month1,day1);
	StringToDate(T2,year2,month2,day2);
	StringToDate(T3,year3,month3,day3);
	int w0 = WeekDate(year0,month0,day0);
	int w1 = WeekDate(year1,month1,day1);
	int w2 = WeekDate(year2,month2,day2);
	int w3 = WeekDate(year3,month3,day3);
	int Total_days = DaysBetween2Date(T0,T1);
	int Total_days_T2 = DaysBetween2Date(T0,T2);
	int Total_days_T3 = DaysBetween2Date(T0,T3);
	
	vector<int> Weekend;
	int F = (6-w0)*m;
	int Lim;
	if(w0 < w1 ){
		 Lim =  (Total_days/7);
	}
	else{
		 Lim = (Total_days/7 +1 );
	}
	
	for(int i=0;i< Lim;i++){
		
		Weekend.push_back(F+1);
		Weekend.push_back(F+2);
		Weekend.push_back(F+3);
		Weekend.push_back(F+4);
		F += 7*2;
	}
	/*for(vector<int>::iterator it=Weekend.begin();it!=Weekend.end();it++){
		
		cout<<"Sephiria The weekend:"<<*it<<endl;
	}*/
	int Trade_days;	
	if(w0 < w1 ){
		Trade_days =Total_days - (Total_days/7)*2 + 1;
	}
	else{
		Trade_days =Total_days - (Total_days/7 +1 )*2 + 1;
	}
	double First,End;
	if(w0 < w2 ){
		First =Total_days_T2 - (Total_days_T2/7)*2;
	}
	else{
		First =Total_days_T2 - (Total_days_T2/7 + 1 )*2;
	}
	if(w0 < w3 ){
		End =Total_days_T3 - (Total_days_T3/7)*2 + 1;
	}
	else{
		End =Total_days_T3 - (Total_days_T3/7 + 1 )*2 + 1;
	}
	int holiday;
	if(w0 < w2 ){
		holiday =(Total_days/7)*2;
	}
	else{
		holiday =(Total_days/7 + 1 )*2;
	}
	
	
	Total_days += 1;
	Total_days_T3 += 1;
	double t2 =Total_days/360.0;
	double n2 = Total_days*m;
/*	double First,End;
	First =Total_days_T2;
	 End =Total_days_T3;*/
	//cout<<"Total_days:"<<End<<endl;
	double t =Trade_days/261.0;
	int n = Trade_days*m;
	/*double t =44.0/261.0;
	double n = 44.0*m;*/
	//	cout<<"GDFG:"<<n<<endl;
	First *=m;
	End*=m;
//cout<< "t"<<setprecision(6)<<t<<endl;
	R = exp(r*(t2/n2));
	double Rinv=1/R;
	
	u = exp(s*sqrt(t/n));
	d= 1/u;
	double uu=u*u;
	double p_up = prob(R,u,d);
	double p_down=1-p_up;
	vector<double> prices(n+1);
	prices[0]=S*pow(d,n)*pow(R,holiday*m);

	for(int i=1;i<=n;++i) prices[i] = uu*prices[i-1];
	vector<double> put_values(n+1);
	//int check = 0;
	//int k =0;
	//int ite=0;
	for(int i=0;i<=n;++i) put_values[i] = max(0.0,(X-prices[i]));
	for(int step=n-1; step>=0; --step){
		/*
		for(int j=Weekend.size()-1;j>=0;j--){
				if (step+1 == Weekend[j]){
					check = 1;
					k=j;
					break; }
			}*/ /*
		if (check == 1){
			for(int i=0;i<=step;++i){
				call_values[i] = call_values[i]*Rinv;
			//call_values[i] = (0.5*call_values[i+1]+0.5*call_values[i]);//Rinv;
	//		prices[i] = d*prices[i+1];
	//	if(step > First&&step <= End)	call_values[i] = max(call_values[i],prices[i]-X);
			}
			check = 3;
		}*/
		/*else if(check == 2){
			for(int i=0;i<=step;++i){
			call_values[i] = (p_up*call_values[i+3]+p_down*call_values[i])*Rinv;	
		//call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
			}
			check = 3;
		}*/
		//else if(check == 3){
		for(int i=0;i<=step;++i){
			put_values[i] = (p_up*put_values[i+1]+p_down*put_values[i])*Rinv;
			prices[i] = d*prices[i+1];
		if(step > First&&step <= End) put_values[i] = max(put_values[i],X-prices[i]);
			if(((n-w1*m)-step) % 10 == 0){
				put_values[i] = put_values[i]*pow(Rinv,2*m);
				 prices[i] = (prices[i])*pow(Rinv,2*m);
				//if(step > First&&step <= End) put_values[i] = max(put_values[i],X-prices[i]);
		}
	}
	}
		return put_values[0];
}

int main(){
	double S,X,s,t,r;
	int m;
	string T0,T1,T2,T3;
	/*cout<<"stock price :";
	cin>>S;
	cout<<"strike price :";
	cin>>X;
	cout<<"annual volatility in percentage :";
	cin>>s;
	cout<<"starting date :";
	cin>>T0;
	cout<<"maturity date of option :";
	cin>>T1;
	cout<<"the number of periods per day for the tree :";
	cin>>m;
	cout<<"interest rate in percentage :";
	cin>>r;
	cout<<"starting date when option can be exercised :";
	cin>>T2;
	cout<<"final date when option can be exercised :";
	cin>>T3;*/ 
	S=100;X=95;s=35;
	T0="2015-04-21";
	T1="2015-06-19";
	m=2;r=10;
	T2="2015-05-21";
	T3="2015-06-19";
	s=s/100;
	r=r/100;
	
	cout<<"The Bermuda options call price :"<<callprice_Bermu(S,X,s,T0,T1,m,r,T2,T3)<<endl;
	cout<<"The Bermuda options put price :"<<putprice_Bermu(S,X,s,T0,T1,m,r,T2,T3)<<endl;
	system("pause");
	return 0;
	
}

