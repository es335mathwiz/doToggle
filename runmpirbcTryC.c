

/*Mathematica Creation Date{2017, 3, 3, 17, 48, 5.631265}*/
/*rbc example model*/
#include "runItExternalDefs.h"
#include "distStochSims.h"
main(int argc, char * argv[])
{
#include "runItInvariantLocalDefs.h"
#include "runrbcTryCLocalDefs.h"
#include "runItInvariantMpiDefs.h"
printf("$Id: mpiRunIt.mc,v 1.1 2001/06/19 19:49:23 m1gsa00 Exp m1gsa00 $
");



  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Get_processor_name(processorName, &nameLength);

/*obtain dimensions for model*/
rbcExampleModelDimensions(&numberOfEquations,&lags,&leads,
	&numberOfParameters,&numberOfDataValues,&numberOfShocks,&numberExog);

/*allocate space for objects that do not depend on command line switches*/

allocLinearTerminator(numberOfEquations,lags,leads,
numberExog,
maxNumberElements,
&upsilonMatrix,&upsilonMatrixj,&upsilonMatrixi,
&hMat,&hMatj,&hMati,
&hzMat,&hzMatj,&hzMati,
&cstar,&cstarj,&cstari,
&AMqMatrix,&AMqMatrixj,&AMqMatrixi,
& rootr,&rooti,
&AMbMatrix,&AMbMatrixj,&AMbMatrixi,
&phiInvMat,&phiInvMatj,&phiInvMati,
&fmat,&fmatj,&fmati,
&varthetaC,&varthetaCj,&varthetaCi,
&varthetaZstar,&varthetaZstarj,&varthetaZstari
);
/*space for data and shocks*/
allocShocksData(numberOfEquations,numberOfShocks,numberOfDataValues,
	&rbcExampleShockVals,&rbcExampleDataVals,
	&rbcExampleZeroShock);
/*space for qmatrix*/
/*allocAltComputeAsymptoticQ(numberOfEquations,lags,leads,spaMaxNumberElements,
	&AMqMatrix,&AMqMatrixj,&AMqMatrixi,&rootr,&rooti);*/
/*space for bmatrix*/
/*allocAltComputeAsymptoticQ(numberOfEquations,lags,leads,spaMaxNumberElements,
	&AMbMatrix,&AMbMatrixj,&AMbMatrixi,&brootr,&brooti);*/
/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after storage allocations
 totalTime=%f,
	userSystemTime=%f,systemTime=%f
",
	*totalTime,*userSystemTime,*(userSystemTime+1));



/*initialize data and shocks*/
for(i=0;i<numberOfDataValues;i++){rbcExampleData(i,
	rbcExampleDataVals+(i*numberOfEquations));}
for(i=0;i<numberOfShocks;i++){rbcExampleShocks(i,
	rbcExampleShockVals+(i*numberOfEquations));}


processCommandLine(argc,argv,namesArray,numberOfEquations,
paramNamesArray,numberOfParameters,parameters,
	rbcExampleDataVals,numberOfDataValues,
	&pathLength,&replications,&t0,&stochasticPathLength,
intControlParameters,doubleControlParameters,flnm);


/*open output file*/
    outFile=fopen(flnm,"w");



/*allocate space for objects that depend on command line switches*/
allocMa50(numberOfEquations,lags,leads,pathLength,pathLength*MAXELEMENTS,
		  &ma50bdIptru,
		  &ma50bdIptrl,
		  &ma50bdIrnf,
		  &ma50bdFact,
		  &ma50bdIq,
		  &ma50bdJob);
allocMa50(numberOfEquations,lags,leads,1,MAXELEMENTS,
		  &cmpma50bdIptru,
		  &cmpma50bdIptrl,
		  &cmpma50bdIrnf,
		  &cmpma50bdFact,
		  &cmpma50bdIq,
		  &cmpma50bdJob);
/*space for FP and for newton step workspace*/
allocFPNewt(numberOfEquations,lags,leads,
	pathLength,MAXELEMENTS,
&rbcExampleFP,
&rbcExampleIntercept,
&fmats,&fmatsj,&fmatsi,&smats,&smatsj,&smatsi);
/*space for path*/
allocPathNewt(numberOfEquations,lags,leads,
	pathLength,replications,stochasticPathLength,
&rbcExamplePathQ,&rbcExampleZeroPathQ);
/*space for stochSims sucess record*/
allocStochSims(stochasticPathLength,replications,&failedQ);
/*initialize  whole path to data values at t0*/
for(i=0;i<lags+pathLength+leads+stochasticPathLength;i++){
  for(j=0;j<numberOfEquations;j++){
	rbcExampleZeroPathQ[i* numberOfEquations+j]=
	  rbcExampleDataVals[(i+t0)*numberOfEquations+j];
	rbcExamplePathQ[i* numberOfEquations+j]=
	  rbcExampleDataVals[(i+t0)*numberOfEquations+j];
  }}



/*initialize  whole path to data values at t0*/
for(i=0;i<lags+1+leads;i++){
  for(j=0;j<numberOfEquations;j++){
	rbcExampleFP[i* numberOfEquations+j]=
	  rbcExampleDataVals[(i+pathLength+t0)*numberOfEquations+j];
  }}








/*
rbcExamplePeriodicPointGuesser(parameters,1,rbcExampleFP);

printf("initiating FP solution computation
");
FPnewt(&numberOfEquations,&lags,&leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,
failedQ,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo);
printf("computed FP solution
");
*/
/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after fixed point computation
 totalTime=%f,userSystemTime=%f,systemTime=%f
",*totalTime,*userSystemTime,*(userSystemTime+1));






printf("generating perm vec
");
allocGenerateDraws(1,stochasticPathLength,replications,
&rbcExamplePermVec);
 generateDraws(1,(stochasticPathLength),replications,numberOfShocks,rbcExamplePermVec);
printf("done generating perm vec
");

/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after generating draws
 totalTime=%f,userSystemTime=%f,
systemTime=%f
",*totalTime,*userSystemTime,*(userSystemTime+1));








/*compute asymptotic Q constraint*/
if(!useIdentityQ){
altComputeAsymptoticQMatrix(
&numberOfEquations,&lags,&leads,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleFP,pathLength,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&spaMaxNumberElements,
AMqMatrix,AMqMatrixj,AMqMatrixi,&auxInit,&qRows,rootr,rooti,
failedQ,0,
 intControlParameters, doubleControlParameters,
 intOutputInfo,  doubleOutputInfo
);
/*
rbcExampleUpsilon(parameters,
upsilonMatrix,upsilonMatrixj,upsilonMatrixi
			  );
linearTerminator(&numberOfEquations,&lags,&leads,
&numberExog,rbcExampleFP,
rbcExample,rbcExampleDerivative,rbcExampleExogH,
parameters,
upsilonMatrix,upsilonMatrixj,upsilonMatrixi,
&spaMaxNumberElements,
hMat,hMatj,hMati,
hzMat,hzMatj,hzMati,
cstar,cstarj,cstari,
AMqMatrix,AMqMatrixj,AMqMatrixi,
&auxInit,&qRows,rootr,rooti,
AMbMatrix,AMbMatrixj,AMbMatrixi,
phiInvMat,phiInvMatj,phiInvMati,
fmat,fmatj,fmati,
varthetaC,varthetaCj,varthetaCi,
varthetaZstar,varthetaZstarj,varthetaZstari,
failedQ
);


*/

fprintf(outFile,"rbcExampleAuxInitQRows={%d,%d}
",auxInit,qRows);
fPrintMathDbl(outFile,(numberOfEquations*(leads+lags)),
      rootr,"rbcExampleRootr");
fPrintMathDbl(outFile,(numberOfEquations*(leads+lags)),
      rooti,"rbcExampleRooti");

/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after computing Q,max elems=%d
 totalTime=%f,userSystemTime=%f,
systemTime=%f
",spaMaxNumberElements,*totalTime,*userSystemTime,*(userSystemTime+1));
}

if(qRows <numberOfEquations*leads){
printf("qmatrix has %d rows, need %d so using identity matrix terminal condition
",qRows,numberOfEquations*leads);
altComputeAsymptoticIMatrix(
&numberOfEquations,&lags,&leads,
AMqMatrix,AMqMatrixj,AMqMatrixi,
failedQ
);
}

/*
spaMaxNumberElements=SPAMAXELEMENTS;
obtainSparseReducedForm(&spaMaxNumberElements,numberOfEquations* leads,
numberOfEquations* (leads+lags),
AMqMatrix,AMqMatrixj,AMqMatrixi,
AMbMatrix,AMbMatrixj,AMbMatrixi);*/
/*time used so far?*/
/* *totalTime=dtime_(userSystemTime);
printf("after computing B
 totalTime=%f,userSystemTime=%f,
systemTime=%f
",*totalTime,*userSystemTime,*(userSystemTime+1));
*/








/*use reduced form to compute rest of path*/
/*for(i=0;i< pathLength+leads+stochasticPathLength;i++){
	applySparseReducedForm(numberOfEquations,
		numberOfEquations* lags,
		rbcExampleZeroPathQ+(i *numberOfEquations),rbcExampleFP,
		AMbMatrix,AMbMatrixj,AMbMatrixi,
rbcExampleZeroPathQ+((i * numberOfEquations)+(numberOfEquations*lags)));
for(j=0;j<numberOfEquations;j++){
rbcExampleZeroPathQ[((i * numberOfEquations)+(numberOfEquations*lags))+j]=
rbcExampleFP[j+(numberOfEquations*lags)%
	(numberOfEquations*(lags+leads+1))];}
}*/

/*set terminal time and call stochSim*/
tf=t0+stochasticPathLength-1;
printf("saving values for variable in file named %s
",flnm);
fprintf(outFile,"rbcExampleRunParams={%d,%d,%d,%d,%d,%d,%d};
",
    numberOfEquations,lags,leads,
     pathLength,t0,stochasticPathLength,replications);
fPrintMathInt(outFile,replications * (stochasticPathLength),
      rbcExamplePermVec,"rbcExamplePermVec");
fPrintMathDbl(outFile,(numberOfEquations*numberOfDataValues),rbcExampleDataVals,"rbcExampleData");
fPrintMathDbl(outFile,(numberOfEquations*numberOfShocks),rbcExampleShockVals,"rbcExampleShocks");



/*begin biiiiiiiiiig MPI paste*/
  if (myRank == 0)
    {
      startwtime = MPI_Wtime();
      printf("Process #%d: %d replications.
", myRank, replications);
    }

  printf("Process %d has been started on host %s
", myRank, processorName);

  pathQLength = replications * numberOfEquations *
    (lags + leads + pathLength + stochasticPathLength);
  failedQLength = replications;
  buildResultType(rbcExamplePathQ, failedQ, &completedDraw, pathQLength,
		  failedQLength, &resultMessageType);

  if ((myRank == 0) && (replications < numberOfProcesses - 1))
    {
      printf("Process #%d: Too many processes for the number of replications. Killing some 
processes.
",
	     myRank);
      for (i = replications + 1; i < numberOfProcesses; i++)
	{
	  sendHaltMessage(i);
	}
      numberOfProcesses = replications + 1;
    }

  if (myRank == 0)
    {
      printf("Process #0: pid = %d
", getpid());
      for (newDraw = 0; newDraw < numberOfProcesses - 1; newDraw++)
	{
	  destination = newDraw + 1;
	  sendDataMessage(destination, &newDraw);
	}

      for (numberOfCompletedDraws = 0; numberOfCompletedDraws < replications;
	   numberOfCompletedDraws++)
	{
	  printf("Process #0: About to start waiting for message
"); fflush(stdout);
	  MPI_Recv(rbcExamplePathQ, 1, resultMessageType, MPI_ANY_SOURCE,
		   RESULT_MSG_TAG, MPI_COMM_WORLD, &status);
	  printf("Process #0: Received message
"); fflush(stdout);
	  source = status.MPI_SOURCE;
	  if (status.MPI_ERROR != 0)
	    error(myRank, source, status.MPI_TAG, status.MPI_ERROR);
	  printf("Process #0 received results of replication %d from process #%d.
",
		 completedDraw, source); fflush(stdout);

	  /**************
	  writeOutput(flnm, completedDraw, rbcExamplePathQ, failedQ, pathQLength,
		      failedQLength, pathLength, t0, stochasticPathLength,
		      replications, rbcExamplePermVec, rbcExampleDataVals,
		      rbcExampleShockVals, lags, leads, numberOfEquations, numberOfShocks,
		      numberOfDataValues);
		      **************/
	  if (numberOfCompletedDraws < replications - numberOfProcesses + 1)
	    {
	      destination = source;
	      sendDataMessage(destination, &newDraw);
	      newDraw++;
	    }
	  else /*** tell other process that we're done ***/
	    {
	      destination = source;
	      sendHaltMessage(destination);
	    }
	}
      printf("Process #%d: All replications have completed.
", myRank); fflush(stdout);
      endwtime = MPI_Wtime();
      printf("Process #%d: Wall clock time = %f.
", myRank, endwtime - startwtime);
      fflush(stdout);
    }
  else /*** myRank not 0 ***/
    {
      int halt = 0;

      printf("Process #%d: pid = %d
", myRank, getpid());

      while (!halt)
	{
	  source = 0;
	  printf("Process #%d: About to start waiting for message
", myRank); fflush(stdout);
	  MPI_Recv(buffer, BUFFER_SIZE, MPI_PACKED, source, MPI_ANY_TAG, MPI_COMM_WORLD,
		   &status);
	  printf("Process #%d: Received message
", myRank); fflush(stdout);
	  printf("Process #%d: About to check for errors. source = %d, tag = %d, status. MPIERROR = 
%d.
",
		 myRank, source, tag, status.MPI_ERROR); fflush(stdout);
	  tag = status.MPI_TAG;
	  if (status.MPI_ERROR != 0)
	    error(myRank, source, tag, status.MPI_ERROR);
	  fflush(stderr);
	  printf("Process #%d: Completed error check.
", myRank); fflush(stdout);

	  if (tag == HALT_MSG_TAG)
	    halt = 1;
	  else
	    {
	      position = 0;
	      MPI_Unpack(buffer, BUFFER_SIZE, &position, &newDraw, 1, MPI_INT,
			 MPI_COMM_WORLD);

	      sprintf(outFileName, "%s%d", flnm, newDraw);
	      outFile=fopen(outFileName,"w");

	      /*set terminal time and call stochSim*/
	      tf=t0+stochasticPathLength-1;

	      printf("Process #%d: saving values for variable in file named 
%s
",myRank,outFileName);
	      fprintf(outFile,"rbcExampleRunParams={%d,%d,%d,%d,%d,%d,%d};
",
		      numberOfEquations,lags,leads,
		      pathLength,t0,stochasticPathLength,replications);
	      fPrintMathInt(outFile,replications * (stochasticPathLength),
			    rbcExamplePermVec,"rbcExamplePermVec");
	      /************************  these generate too much output for frbus
	      
fPrintMathDbl(outFile,(numberOfEquations*numberOfDataValues),rbcExampleDataVals,"rbcExampleData");
	      
fPrintMathDbl(outFile,(numberOfEquations*numberOfShocks),rbcExampleShockVals,"rbcExampleShocks");
/*end biiiiiiiiiig MPI paste*/

*ma50bdJob=1;
pathNewt(&numberOfEquations,&lags,&leads,&pathLength,
rbcExample,rbcExampleDerivative,parameters,
rbcExampleZeroShock,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
rbcExampleFP,rbcExampleIntercept,rbcExampleFP,
rbcExampleZeroPathQ,
failedQ,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru);
fPrintMathInt(outFile,widthIntOutputInfo,intOutputInfo,
"rbcExampleZeroShockIntOutputInfo");
fPrintMathDbl(outFile,widthIntOutputInfo,doubleOutputInfo,
"rbcExampleZeroShockDoubleOutputInfo");


fPrintMathDbl(outFile,(numberOfEquations*(leads+pathLength+lags)),
      rbcExampleZeroPathQ,"rbcExampleZeroShockResults");
fPrintMathInt(outFile,1,
failedQ,"rbcExampleZeroShocksFailedQ");
/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after using Q matrix
totalTime=%f,userSystemTime=%f,
systemTime=%f
",*totalTime,*userSystemTime,*(userSystemTime+1));
*cmpma50bdJob=1;
diststochSim(&numberOfEquations,&lags,&leads,&pathLength,
rbcExample,rbcExampleDerivative,parameters,
&replications,&t0,&tf,rbcExamplePermVec,
rbcExampleShockVals,&numberOfShocks,
rbcExampleDataVals,&numberOfDataValues,
fmats,fmatsj,fmatsi,
smats,smatsj,smatsi,
&maxNumberElements,AMqMatrix,AMqMatrixj,AMqMatrixi,
rbcExampleZeroPathQ,rbcExampleIntercept,rbcExampleFP,
rbcExamplePathQ,
failedQ,intControlParameters,doubleControlParameters,
intOutputInfo, doubleOutputInfo,
ma50bdJob,
ma50bdIq,
ma50bdFact,
ma50bdIrnf,
ma50bdIptrl,
ma50bdIptru,
cmpma50bdJob,
cmpma50bdIq,
cmpma50bdFact,
cmpma50bdIrnf,
cmpma50bdIptrl,
cmpma50bdIptru,newDraw,outFile
);
fPrintMathDbl(outFile,(replications * numberOfEquations*(stochasticPathLength+lags)),
      rbcExamplePathQ,"rbcExampleResults");
fPrintMathInt(outFile,replications*widthIntOutputInfo,intOutputInfo,
"rbcExampleIntOutputInfo");
fPrintMathDbl(outFile,replications*widthIntOutputInfo,doubleOutputInfo,
"rbcExampleDoubleOutputInfo");
fPrintMathInt(outFile,replications*stochasticPathLength,
failedQ,"rbcExamplefailedQ");
/*time used so far?*/
*totalTime=dtime_(userSystemTime);
printf("after using Q matrix
totalTime=%f,userSystemTime=%f,
systemTime=%f
",*totalTime,*userSystemTime,*(userSystemTime+1));
	      writeOutput(outFile, newDraw, rbcExamplePathQ, failedQ, pathQLength,
			  failedQLength, pathLength, t0, stochasticPathLength,
			  replications, rbcExamplePermVec, rbcExampleDataVals,
			  rbcExampleShockVals, lags, leads, numberOfEquations, numberOfShocks,
			  numberOfDataValues);

	      printf("Process #%d: Finished stochSim().
", myRank); fflush(stdout);

	      destination = 0;
	      tag = RESULT_MSG_TAG;
	      completedDraw = newDraw;
	      printf("Process #%d sending replication %d results to process #0.
",
		     myRank, completedDraw);
	      MPI_Send(rbcExamplePathQ, 1, resultMessageType, destination,
		       tag, MPI_COMM_WORLD);
	      fclose(outFile);
	    }
}}

  MPI_Finalize();






cfreeLinearTerminator(
&upsilonMatrix,&upsilonMatrixj,&upsilonMatrixi,
&hMat,&hMatj,&hMati,
&hzMat,&hzMatj,&hzMati,
&cstar,&cstarj,&cstari,
&AMqMatrix,&AMqMatrixj,&AMqMatrixi,
& rootr,&rooti,
&AMbMatrix,&AMbMatrixj,&AMbMatrixi,
&phiInvMat,&phiInvMatj,&phiInvMati,
&fmat,&fmatj,&fmati,
&varthetaC,&varthetaCj,&varthetaCi,
&varthetaZstar,&varthetaZstarj,&varthetaZstari
);

cfreeGenerateDraws(&rbcExamplePermVec);
cfreeShocksData(&rbcExampleShockVals,&rbcExampleDataVals,
	&rbcExampleZeroShock);
/*cfreeAltComputeAsymptoticQ(
&AMqMatrix,&AMqMatrixj,&AMqMatrixi,&rootr,&rooti);
cfreeAltComputeAsymptoticQ(
&AMbMatrix,&AMbMatrixj,&AMbMatrixi,&brootr,&brooti);*/
cfreePathNewt(&rbcExamplePathQ);
cfreePathNewt(&rbcExampleZeroPathQ);
cfreeFPNewt(lags,pathLength,
&rbcExampleFP,
&rbcExampleIntercept,
&fmats,&fmatsj,&fmatsi,&smats,&smatsj,&smatsi);
cfreeStochSims(&failedQ);
cfreeMa50(&ma50bdIptru,
		  &ma50bdIptrl,
		  &ma50bdIrnf,
		  &ma50bdFact,
		  &ma50bdIq,
		  &ma50bdJob);
cfreeMa50(&cmpma50bdIptru,
		  &cmpma50bdIptrl,
		  &cmpma50bdIrnf,
		  &cmpma50bdFact,
		  &cmpma50bdIq,
		  &cmpma50bdJob);


return(0);

}

#define SHOCKS 4
void modData(int numberOfEquations,int numberDataValues,double * dataVals,
			 int vbl,int t0,int tf,double val1,double val2)
{
  int t;
  for(t=t0;t<=tf&&t<numberDataValues;t++){
dataVals[t*numberOfEquations+vbl]=dataVals[t*numberOfEquations+vbl]+
  (t-t0)*val2/(tf-t0) + (tf-t)*val1/(tf-t0);
  }
}
void modDataAbs(int numberOfEquations,int numberDataValues,double * dataVals,
			 int vbl,int t0,int tf,double val1,double val2)
{
  int t;
  for(t=t0;t<=tf&&t<numberDataValues;t++){
dataVals[t*numberOfEquations+vbl]=(t-t0)*val2/(tf-t0) + (tf-t)*val1/(tf-t0);
  }
}
 
 
#include "runItOther.h"

