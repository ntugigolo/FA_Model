#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <limits>

using namespace std;
double Price(double y,vector<double> c,double w){
	int Size=c.size();
	double price=0;
	for (int i=0;i<Size;i++){
		price+=c[i]/pow((1+y),w+i);  
	}
	return price;
}
void Duration(double y,vector<double> c,double w,double p){
	int Size=c.size();
	double duration=0;
	for (int i=0;i<Size;i++){
		duration+=((w+i)*c[i]/pow(1+y,w+i+1))/p; 
		}
		cout<<"The modified duration is:"<<duration<<endl;
	}
void Convexity(double y,vector<double> c,double w,double p){
	int Size=c.size();
	double convexity=0;
	for (int i=0;i<Size;i++){
		convexity+=((w+i)*(w+i+1)*c[i]/pow(1+y,w+i+2))/p;
			}
		cout<<"The convexity is:"<<convexity<<endl;
	 }
void Check(int c,string a){
	 string::size_type pos,pos_check,pos_check_cash;
	 string check="0123456789.",check_cash="01234567890.,";
	 double W = atof(a.c_str());
	 int j=0;
     int size=a.size();
	 int size_check=check.size(),size_cash=check_cash.size();
		while(j<size){
			pos_check=check.find(a[j],0);
			pos_check_cash=check_cash.find(a[j],0);
			if(c!=3){
				if(pos_check>size_check){
					switch(c){
						case 1:
							cout<<"先生/女士，您的Interest Rate資料輸入格式錯誤喔!"<<endl;
							system("pause");
							exit(1);
						case 2:
							cout<<"先生/女士，您的w資料輸入格式錯誤喔!"<<endl;
							system("pause");
							exit(1);
							}
				}
			}
			if(c==3) {
				if(pos_check_cash>size_cash){
					cout<<"先生/女士，您的Cash flow資料輸入格式錯誤喔!"<<endl;
					system("pause");
					exit(1);
					}
				}
			else{ 
				if(W>1&&c==2){
					cout<<"W 的範圍應在0~1之間"<<endl;
					system("pause");
					exit(1);
						}
					}	
		j++;}
		}
vector<string> split(string str,string seperate){
     string::size_type pos,pos_cash;
     vector<string> result;
     str+=seperate;
     int size=str.size();
	 for(int i=0;i<size;i++){
         pos=str.find(seperate,i);
		 if(pos==i){
			cout<<"先生/女士，您的Cash flow資料輸入格式錯誤喔!"<<endl;
			system("pause");
			exit(1);
				}
         if(pos<size)
         {
             string s=str.substr(i,pos-i);
             result.push_back(s);
             i=pos+seperate.size()-1;
         }
     }
     return result;
 }
int main(){
	string c,y,w,d = ",";
	vector<double> Cashflow;
	cout<<"Interst Rate:";
	cin>>y;
	Check(1,y);
	double Y = atof(y.c_str());
	cin.ignore (std::numeric_limits<std::streamsize>::max(), '\n'); 
	cout<<"請輸入Cash Flow(以逗號為間隔 ex:1,2,3,101):";
	getline(cin,c);
	Check(3,c);
	vector<string> Result = split(c,d);
	cout<<"W:";
	cin>>w;
	Check(2,w);
	double W = atof(w.c_str());
	for (int i=0;i<Result.size();i++){
	double cashflow = atof(Result[i].c_str());
	Cashflow.push_back(cashflow);
	}
	double p = Price(Y,Cashflow,W);
	Duration(Y,Cashflow,W,p);
	Convexity(Y,Cashflow,W,p);
	system("pause");
	return 0;
}

