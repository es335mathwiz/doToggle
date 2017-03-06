
/*Mathematica Creation Date {2017, 3, 3, 17, 48, 5.631265}*/
/*rbc example model*/
#include "/msu/home/m1gsa00/git/stackStochSims/lagLead.h"
#include <math.h>
#define aDummy(t)     (stateVector[(t-(-1))*4+0])
#define cc(t)     (stateVector[(t-(-1))*4+1])
#define kk(t)     (stateVector[(t-(-1))*4+2])
#define theta(t)     (stateVector[(t-(-1))*4+3])
#define linPt$aDummy(t)     (linearizationPoint[(t-(-1))*4+0])
#define linPt$cc(t)     (linearizationPoint[(t-(-1))*4+1])
#define linPt$kk(t)     (linearizationPoint[(t-(-1))*4+2])
#define linPt$theta(t)     (linearizationPoint[(t-(-1))*4+3])

#define modelShock(n) (shockVec[n])


void rbcExample(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat)
{
double homotopyAlpha[1]={1.0};double linearizationPoint[1]={0.0};
rbcExampleHomotopy(stateVector,parameters,shockVec,aMat,jaMat,iaMat,homotopyAlpha,linearizationPoint);
}


void rbcExampleHomotopy(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat,double * homotopyAlpha,double * linearizationPoint
)
{
int i;
double bMat[4];
int ibMat[4+1];
int jbMat[4];
if(*homotopyAlpha>=1.0) {
double okay1;
double okay5;
double okay7;


okay1=cc(0);

okay5=kk(0);

okay7=theta(0);

aMat[0]=1/okay1-(0.34199999999999997*okay7*pow(okay5,-0.64))/cc(tPaOne)\
;

aMat[1]=okay1+okay5-okay7*pow(kk(tMaOne),0.36);

aMat[2]=okay7-pow(theta(tMaOne),0.95);

aMat[3]=aDummy(0);

/*for(i=0;i<4-0;i++){aMat[i]=aMat[i]+shockVec[i];};*/
iaMat[0]=1.;

iaMat[1]=2.;

iaMat[2]=3.;

iaMat[3]=4.;

iaMat[4]=5.;

jaMat[0]=1.;

jaMat[1]=1.;

jaMat[2]=1.;

jaMat[3]=1.;

} else {
double okay1;
double okay10;
double okay11;
double okay15;
double okay2;
double okay21;
double okay24;
double okay25;
double okay26;
double okay30;
double okay35;
double okay39;
double okay8;
double okay9;


okay2=linPt$cc(0);

okay11=linPt$kk(0);

okay8=linPt$cc(tPaOne);

okay15=linPt$theta(0);

okay9=1/okay8;

okay21=pow(okay11,-0.64);

okay1=cc(0);

okay10=kk(0);

okay30=linPt$kk(tMaOne);

okay35=pow(okay30,0.36);

okay24=-okay15;

okay25=theta(0);

okay26=okay24+okay25;

okay39=linPt$theta(tMaOne);

aMat[0]=1/okay2-0.34199999999999997*okay15*okay21*okay9-0.3419999999999\
9997*okay21*okay26*okay9+0.21888*(okay10-okay11)*okay15*okay9*pow(okay\
11,-1.6400000000000001)-(okay1-okay2)*pow(okay2,-2.)+0.341999999999999\
97*okay15*okay21*(-okay8+cc(tPaOne))*pow(okay8,-2.);

aMat[1]=okay1+okay10-okay15*okay35-okay26*okay35-0.36*okay15*(-okay30+k\
k(tMaOne))*pow(okay30,-0.64);

aMat[2]=okay25-pow(okay39,0.95)-0.95*pow(okay39,-0.050000000000000044)*\
(-okay39+theta(tMaOne));

aMat[3]=aDummy(0);

/*for(i=0;i<4-0;i++){aMat[i]=aMat[i]+shockVec[i];};*/
jaMat[0]=1.;

jaMat[1]=2.;

jaMat[2]=3.;

jaMat[3]=4.;

jaMat[4]=5.;

iaMat[0]=1.;

iaMat[1]=1.;

iaMat[2]=1.;

iaMat[3]=1.;

if(*homotopyAlpha>0.0) {
double okay1;
double okay10;
double okay11;
double okay12;
double okay13;
double okay17;
double okay19;
double okay23;
double okay28;
double okay3;
double okay30;
double okay31;
double okay34;
double okay35;
double okay40;
double okay46;
double okay48;


okay1=cc(0);

okay3=linPt$cc(0);

okay13=linPt$kk(0);

okay10=linPt$cc(tPaOne);

okay17=linPt$theta(0);

okay11=1/okay10;

okay23=pow(okay13,-0.64);

okay19=cc(tPaOne);

okay12=kk(0);

okay28=theta(0);

okay35=linPt$kk(tMaOne);

okay34=kk(tMaOne);

okay40=pow(okay35,0.36);

okay30=-okay17;

okay31=okay28+okay30;

okay46=linPt$theta(tMaOne);

okay48=theta(tMaOne);

aMat[0]=1/okay1+0.34199999999999997*okay11*okay17*okay23-1/okay3+0.3419\
9999999999997*okay11*okay23*okay31-0.34199999999999997*okay17*(-okay10\
+okay19)*okay23*pow(okay10,-2.)-(0.34199999999999997*okay28*pow(okay12\
,-0.64))/okay19-0.21888*okay11*(okay12-okay13)*okay17*pow(okay13,-1.64\
00000000000001)+(okay1-okay3)*pow(okay3,-2.);

aMat[1]=okay17*okay40+okay31*okay40-okay28*pow(okay34,0.36)+0.36*okay17\
*(okay34-okay35)*pow(okay35,-0.64);

aMat[2]=0.95*(-okay46+okay48)*pow(okay46,-0.050000000000000044)+pow(oka\
y46,0.95)-pow(okay48,0.95);

aMat[3]=0.;

jaMat[0]=1.;

jaMat[1]=2.;

jaMat[2]=3.;

jaMat[3]=4.;

jaMat[4]=5.;

iaMat[0]=1.;

iaMat[1]=1.;

iaMat[2]=1.;

iaMat[3]=1.;

for(i=0;i<4;i++){aMat[i]=aMat[i]+(*homotopyAlpha*bMat[i]);};
}
}
}
