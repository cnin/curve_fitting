#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
using namespace std;

//最小二乘拟合相关函数定义
double sum(vector<double> Vnum, int n);
double MutilSum(vector<double> x, vector<double> y, int dataSize);
double RelatePow(vector<double> x, int n, int ex);
double RelateMutiXY(vector<double> x, vector<double> y, int dataSize, int ex);
void EMatrix(vector<double> x, vector<double> y, int dataSize, int ex, double coefficient[]);
void CalEquation(int exp, double coefficient[]);
double F(double c[],int l,int m);
double Em[6][4];
void CalcCoefA(vector<double> x, vector<double> y, int dataSize, double &coefA);

//主函数，这里将数据拟合成二次曲线
int main(int argc, char* argv[])
{
    FILE* fp=nullptr;
	fp=fopen("data.txt","rb+");
	double arry1[5],arry2[5],arry3[100],arry4[100];
	double coefficient[5];
    memset(coefficient,0,sizeof(double)*5);
    vector<double> x,y;

    for(int i=0;i<5;i++){
	    fscanf(fp,"%lf %lf",&arry1[i],&arry2[i]);
        x.push_back(arry1[i]);
        y.push_back(arry2[i]);
    }
    

    EMatrix(x,y,x.size(),3,coefficient);
    printf("Quadratic Fitting: y = %lf + %lf x + %lf x^2 \n",coefficient[0],coefficient[1],coefficient[2]);
	double coefA;
	CalcCoefA(x,y,x.size(),coefA);
    printf("Linear Fitting: y = %3.10le x \n",coefA);

    FILE* fp_write=nullptr;
    fp_write=fopen("data_result.txt","w+");
	for(int i=0;i<100;i++){
	    arry3[i] = coefficient[0]+coefficient[1]*i/20+coefficient[2]*i*i/20/20;
        arry4[i] = coefA*i/20;
        fprintf(fp_write,"%lf %lf %lf\n",double(i)/20,arry3[i],arry4[i]);
    }

    fclose(fp);
	fclose(fp_write);
    return 0;
}
//一次项拟合
void CalcCoefA(vector<double> x, vector<double> y, int n, double &coefA)
{
    coefA=MutilSum(x,y,x.size())/RelatePow(x,x.size(),2);
}
//累加
double sum(vector<double> Vnum, int n)
{
    double dsum=0;
    for (int i=0; i<n; i++)
    {
        dsum+=Vnum[i];
    }
    return dsum;
}
//乘积和
double MutilSum(vector<double> x, vector<double> y, int dataSize)
{
    double dMultiSum=0;
    for (int i=0; i<dataSize; i++)
    {
        dMultiSum+=x[i]*y[i];
    }
    return dMultiSum;
}
//ex次方和
double RelatePow(vector<double> x, int dataSize, int ex)
{
    double ReSum=0;
    for (int i=0; i<=dataSize; i++)
    {
        ReSum+=pow(x[i],ex);
    }
    return ReSum;
}
//x的ex次方与y的乘积的累加
double RelateMutiXY(vector<double> x, vector<double> y, int dataSize, int ex)
{
    double dReMultiSum=0;
    for (int i=0; i<dataSize; i++)
    {
        dReMultiSum+=pow(x[i],ex)*y[i];
    }
    return dReMultiSum;
}
//计算方程组的增广矩阵
void EMatrix(vector<double> x, vector<double> y, int dataSize, int ex, double coefficient[])
{
    // for (int i=1; i<=ex; i++)
    // {
    //     for (int j=1; j<=ex; j++)
    //     {
    //         Em[i][j]=RelatePow(x,dataSize,i+j-2);
    //     }
    //     Em[i][ex+1]=RelateMutiXY(x,y,dataSize,i-1);


    // }
   
    // Em[1][1]=dataSize;
// //	if(ex<3){
// //	    for(int i=1;i<=ex;i++){
// //		    for(int j=1;j<=ex;j++){
// //			if(i!=2&&j!=2)Em[i][j]=0;}
// //		}
// //
// //	} 


    for (int i=0; i<ex; i++)
    {
        for (int j=0; j<ex; j++)
        {
            Em[i][j]=RelatePow(x,dataSize,i+j);
        }
        Em[i][ex]=RelateMutiXY(x,y,dataSize,i);
    }
   
    Em[0][0]=dataSize;
// 	if(ex<3){
// 	    for(int i=1;i<=ex;i++){
// 		    for(int j=1;j<=ex;j++){
// 			if(i!=2&&j!=2)Em[i][j]=0;}
// 		}

// 	} 

for (int i=0;i<ex;i++){
    for (int j=0;j<=ex;j++){
        printf("(%d %d), %lf\n",i,j,Em[i][j]);
    }

}

    CalEquation(ex,coefficient);
}
//求解方程
void CalEquation(int exp, double coefficient[])
{
    for(int k=0;k<exp-1;k++) //消元过程
    {
        for(int i=k+1;i<=exp;i++)
        {
            double p1=0;

            if(Em[k][k]!=0)
                p1=Em[i][k]/Em[k][k];

            for(int j=k;j<=exp;j++)
                Em[i][j]=Em[i][j]-Em[k][j]*p1;
        }
    }
    for (int i=0;i<exp;i++){
        for (int j=0;j<=exp;j++){
           printf("(%d %d), %lf\n",i,j,Em[i][j]);
        }
    }
    coefficient[exp-1]=Em[exp-1][exp]/Em[exp-1][exp-1];
    printf("coef[%d]=%lf\n",exp-1,coefficient[exp-1]);
    for(int l=exp-2;l>=0;l--){   //回代求解
        coefficient[l]=(Em[l][exp]-F(coefficient,l+1,exp))/Em[l][l];
     //   printf("coef[%d]=%lf\n",l,coefficient[l]);
    }
}
//供CalEquation函数调用
double F(double c[],int l,int m)
{
    double sum=0;
    for(int i=l;i<=m;i++)
        sum+=Em[l-1][i]*c[i];
        printf("%lf\n",sum);
    return sum;
}
