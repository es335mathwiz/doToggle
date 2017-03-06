
/*Mathematica Creation Date{2016, 5, 15, 18, 4, 40.611484}*/
/*rbc example model*/
#include "/Users/garyanderson/git/stackStochSims/lagLead.h"
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

  







void rbcExampleDerivative(double *stateVector,double *parameters,
double * shockVec,
double * aMat,int * jaMat,int *iaMat,double * homotopyAlpha,double * linearizationPoint
)
{int i;
double bMat[11];
int ibMat[5];
int jbMat[11];
double cMat[11];
int icMat[5];
int jcMat[11];
int aOne=1;int ierr;int maxNumberHElements;
int hrows=4;
int hcols=12;
double okay[50000];
if(*homotopyAlpha>=1.0) {




aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;

iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
} else {


aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;

iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
/*initialize cMat to zero sparse matrix*/
for(i=0;i<=5;i++){icMat[i]=1;}
for(i=0;i<11;i++){cMat[i]=0;};
if(*homotopyAlpha>0.0) {


aMat[0]=1.;
aMat[1]=1.;
aMat[2]=1.;
aMat[3]=1.;

iaMat[0]=1.;
iaMat[1]=1.;
iaMat[2]=1.;
iaMat[3]=1.;
jaMat[0]=1.;
jaMat[1]=2.;
jaMat[2]=3.;
jaMat[3]=4.;
jaMat[4]=5.;
for(i=0;i<11;i++){cMat[i]=cMat[i]*(*homotopyAlpha);};
}
maxNumberHElements=11;
aplb_(&hrows,&hcols,&aOne,bMat,jbMat,ibMat,cMat,jcMat,icMat,
aMat,jaMat,iaMat,&maxNumberHElements,okay,&ierr);
}
}
