/*Mathematica Creation Date{2017, 3, 3, 17, 48, 5.631265}*/
/*rbc example model*/
#include <stdlib.h>
#include<stdio.h>

void rbcExampleData(int t,double * vectorOfVals);
void rbcExampleShocks(int t,double * vectorOfVals);
void rbcExamplePeriodicPointGuesser(double * parameters,int period,double*);
/*#include "useSparseAMA.h"*/
#include "stackC.h"
#include "stochProto.h"
/*modelDimensions call determines these*/
FILE * outFile;

#define PATHLENGTH 1000
#define rbcNLAGS 1
#define rbcNLEADS 1
#define rbcNEQS 4
#define SHOCKS 30
#define DATA 50
#define PATHLENGTH 1000

int  numberOfEquations[1]={rbcNEQS};
int  lags[1]={rbcNLAGS};
int  leads[1]={rbcNLEADS};
int pathLength[1]={PATHLENGTH};
int stochasticPathLength[1]={PATHLENGTH};
FILE * outFile;
 /* char outFileName[250];*/
static char flnm[50] = "stochOut.m";
/*a counter*/
int i;/*int j;*/
/*modelDimensions call determines these*/
int t0[1]={0};
int tf[1]={0};
int replications[1]={1};
double totalTime[1];
double userSystemTime[2];
/*int shockIndex[1];*/
int  numberOfParameters;
/*int  numberOfDataValues;*/
/*int  numberOfShocks;*/
/*int  numberExog;*/



int numberOfData[1]={DATA};
int numberOfShocks[1]={SHOCKS};

char * namesArray[] =  {"aDummy", "cc", "kk", "theta"};
char * paramNamesArray[] = {};
int numberOfParameters=0;
int * parameters[]={};
int numDATA=500;
int numSHOCKS=500;
double * theData;


void cfree(void * ptr){free(ptr);}
int main(int argc, char * argv[])
{
#include "runrbcTryCLocalDefs.h"
printf(" runIt.mc, 2016 m1gsa00 \n");

rbcExampleDataVals=(double *)calloc(*numberOfEquations*numDATA,sizeof(double));
for(i=0;i<numDATA;i++){rbcExampleData(i,rbcExampleDataVals+(i**numberOfEquations));}

rbcExampleShockVals=(double *)calloc(*numberOfEquations*numSHOCKS,sizeof(double));
for(i=0;i<numSHOCKS;i++){rbcExampleShocks(i,rbcExampleShockVals+(i**numberOfEquations));}



rbcExampleFP=(double *)calloc(rbcNEQS*(rbcNLEADS+rbcNLAGS+1),sizeof(double));

AMqMatrix=(double *)
   calloc(MAXELEMENTS,sizeof(double));
AMqMatrixj=(int *)
   calloc(MAXELEMENTS,sizeof(int));
AMqMatrixi=(int *)
   calloc((rbcNEQS*(rbcNLEADS+rbcNLAGS)),
        sizeof(int));


processCommandLine(argc,argv,namesArray,*numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	rbcExampleDataVals,numDATA,numSHOCKS,
	pathLength,replications,t0,stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);

    outFile=fopen(flnm,"w");
int * failedQ;

failedQ=(int *)calloc(*replications,sizeof(int));
for(i=0;i<*replications;i++)failedQ[i]=0;
*tf=(*t0)+ *stochasticPathLength-1;
 
 


rbcExamplePermVec=(int *)calloc(
     (*stochasticPathLength)*(*replications),sizeof(int));
rbcExamplePathQ=(double *)calloc(
    *replications*
    rbcNEQS*(rbcNLAGS+rbcNLEADS+(*pathLength)+*stochasticPathLength),
sizeof(double));
double ** ptrToPtrToDouble = NULL;

fmats =(double **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsj =(int **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
fmatsi =(int **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
smats =(double **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
smatsj =(int **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
smatsi =(int **)calloc((*pathLength)+rbcNLAGS+1,sizeof(ptrToPtrToDouble));
for(i=0;i<(*pathLength)+rbcNLAGS+1;i++){
fmats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
fmatsj[i] =(int *)calloc(MAXELEMENTS,sizeof(int));
fmatsi[i] =(int *)calloc(
     rbcNEQS*(rbcNLAGS+rbcNLEADS)+1,sizeof(int));
smats[i] =(double *)calloc(MAXELEMENTS,sizeof(double));
smatsj[i] =(int *)calloc(MAXELEMENTS,sizeof(int));
smatsi[i] =(int *)calloc(
     rbcNEQS*(rbcNLAGS+rbcNLEADS)+1,sizeof(int));
}




rbcExamplePeriodicPointGuesser(parameters,1,rbcExampleFP);

FPnewt(numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,
failedQ);

if(failedQ[0])
{  printf("problems with finding fixedpoint");return(1);}
else {printf("computed fixed point");}






printf("generating perm vec");
 generateDraws(1,(*stochasticPathLength),(*replications),numSHOCKS,rbcExamplePermVec);
printf("done generating perm vec");






altComputeAsymptoticQMatrix(
numberOfEquations,lags,leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,
failedQ
);



if(failedQ[0])
{  printf("problems computing  Q matrix");return(1);}
else {printf("computed Q matrix");}

printf("computed Q matrix");
for(i=0;i< *pathLength;i++){
rbcExamplePeriodicPointGuesser(parameters,1,
rbcExamplePathQ+(i *rbcNEQS));}


/*int  dtime(double * userSystemTime);*/
/**totalTime=dtime(userSystemTime);*/
printf("after computing Q matrix totalTime=%f,userSystemTime=%f,systemTime=%f",*totalTime,*userSystemTime,*(userSystemTime+1));
printf("using q matrix");


stochSim(numberOfEquations,lags,leads,pathLength,
rbcExample,rbcExampleDerivative,parameters,
replications,t0,tf,rbcExamplePermVec,
rbcExampleShockVals,numberOfShocks,
rbcExampleDataVals,numberOfData,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
rbcExampleFP,
rbcExamplePathQ,
failedQ);




/**totalTime=dtime(userSystemTime);*/
clock_t start,end;
double cpu_time_used;
start=clock();
end=clock();
cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
printf("start=%lu,end=%lu,used=%e\n",start,end,cpu_time_used);
#include <unistd.h> 
#include <stdint.h> 
#include <sys/times.h>
#include <sys/types.h>
struct tms t;
printf("tick frequency is: %lu\n", sysconf(_SC_CLK_TCK)); 
if (times(&t) < 0) 
{ perror("times"); /* error - print a message and exit */ exit(1); } 
/* print the results */ 
printf("user time: %ju ticks\n", (uintmax_t) t.tms_utime); 
printf("system time: %ju ticks\n", (uintmax_t) t.tms_stime); 
printf("chidren - user time: %ju ticks\n", (uintmax_t) t.tms_cutime); 
printf("chidren - system time: %ju ticks\n", (uintmax_t) t.tms_cstime); 


printf("after using Q matrix\ntotalTime=%f,userSystemTime=%f,systemTime=%f\n",
*totalTime,*userSystemTime,*(userSystemTime+1));








printf("saving values for variable in file named %s\n",flnm);
fprintf(outFile,"rbcExampleRunParams={%d,%d,%d,%d,%d,%d,%d};\n",
    rbcNEQS,rbcNLAGS,rbcNLEADS,
     *pathLength,*t0,*stochasticPathLength,*replications);
fPrintMathInt(outFile,*replications,failedQ,"rbcExampleFailedQ");
fPrintMathInt(outFile,*replications * (*stochasticPathLength),
      rbcExamplePermVec,"rbcExamplePermVec");
fPrintMathDbl(outFile,(*replications * rbcNEQS*(*stochasticPathLength+rbcNLAGS)),
      rbcExamplePathQ,"rbcExampleResults");
fPrintMathDbl(outFile,(rbcNEQS*(DATA)),rbcExampleDataVals,"rbcExampleData");
fPrintMathDbl(outFile,(rbcNEQS*(SHOCKS)),rbcExampleShockVals,"rbcExampleShocks");
     fclose(outFile);








cfree(failedQ);
cfree(AMqMatrix);
cfree(AMqMatrixj);
cfree(AMqMatrixi);
cfree(rbcExampleShockVals);
cfree(rbcExampleDataVals);
cfree(rbcExamplePermVec);
cfree(rbcExampleFP);
cfree(rbcExamplePathQ);
for(i=0;i<(*pathLength)+rbcNLAGS+1;i++){
cfree(smats[i]);
cfree(smatsj[i]);
cfree(smatsi[i]);
cfree(fmats[i]);
cfree(fmatsj[i]);
cfree(fmatsi[i]);
}
cfree(fmats);
cfree(fmatsj);
cfree(fmatsi);
cfree(smats);
cfree(smatsj);
cfree(smatsi);


 
return(0);
}



