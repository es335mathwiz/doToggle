/*Mathematica Creation Date{2016, 5, 15, 18, 4, 40.665032}*/
/*rbc example model*/
#include<stdlib.h>

#include "/Users/garyanderson/git/stackStochSims/runItExternalDefs.h"



#define PATHLENGTH 1000

int numberOfEquations=4;
char * namesArray[] =  {"aDummy", "cc", "kk", "theta"};
char * paramNamesArray[] = {};;
int numberOfParameters=0;
int * parameters[]={};
int numDATA=500;
int numSHOCKS=500;
double * theData;

void  cfree(void * ptr){free(ptr);}

int main(int argc, char * argv[])
{
#include "/Users/garyanderson/git/stackStochSims/runItInvariantLocalDefs.h"
#include "runrbcTryCLocalDefs.h"
printf("$Id: runIt.mc, 2016 m1gsa00 $\n");

rbcExampleDataVals=(double *)calloc(numberOfEquations*numDATA,sizeof(double));
for(i=0;i<numDATA;i++){rbcExampleData(i,rbcExampleDataVals+(i*numberOfEquations));}

rbcExampleShockVals=(double *)calloc(numberOfEquations*numSHOCKS,sizeof(double));
for(i=0;i<numSHOCKS;i++){rbcExampleShocks(i,rbcExampleShockVals+(i*numberOfEquations));}


processCommandLine(argc,argv,namesArray,numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	rbcExampleDataVals,numDATA,numSHOCKS,
	&pathLength,&replications,&t0,&stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);

/*
rbcExamplePeriodicPointGuesser(parameters,1,rbcExampleFP);

FPnewt(numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
failedQ);

*/
}

#include "/Users/garyanderson/git/stackStochSims/runItOther.h"

/*
printf("generating perm vec\n");
 generateDraws(1,(stochasticPathLength),(*replications),numSHOCKS,julliardPermVec);
printf("done generating perm vec\n");
*/


/*


altComputeAsymptoticQMatrix(
numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,
failedQ
);
*/

/*
if(failedQ[0])
{  printf("problems computing  Q matrix\n");return(1);}
else {printf("computed Q matrix\n");}
*//*
printf("computed Q matrix\n");
for(i=0;i< *pathLength;i++){
rbcExamplePeriodicPointGuesser(parameters,1,
julliardPathQ+(i *julNEQS));}
*totalTime=dtime(userSystemTime);
printf("after computing Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix\n");*/
/*

stochSim(numberOfEquations,lags,leads,pathLength,
rbcExample,rbcExampleDerivative,parameters,
replications,t0,tf,julliardPermVec,
julliardShocks,numberOfShocks,
julliardData,numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
julliardFP,
julliardPathQ,
failedQ);

*/

