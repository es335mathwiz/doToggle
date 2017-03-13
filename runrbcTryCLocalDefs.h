



/*Mathematica Creation Date{2017, 3, 3, 17, 48, 5.631265}*/
#define SPAMAXELEMENTS 15*(4)^2
#define MAXELEMENTS 25*(4)^2

int  maxNumberElements=MAXELEMENTS;
/*int  spaMaxNumberElements=SPAMAXELEMENTS;*/
void rbcExample(double * xvec,double * pvec,double * shock,
double * alhs,
int * jalhs,int * ialhs,int * alphas,double * linPt
);
void rbcExampleDerivative(double * xvec,double * pvec,
double * alhs,
int * jalhs,
int * ialhs);
void rbcExampleExogH(double * pvec,
double * alhs,
int * jalhs,
int * ialhs);
/*model specific names and data*/
char * namesArray[] =  
{"aDummy", "cc", "kk", "theta"};
char * paramNamesArray[] =  
{};
double parameters[]=
{};
/*int rbcExampleexogQ[]=
{0, 0, 0, 0};*/

int * rbcExamplePermVec;
/*double * rbcExampleZeroShock;*/
double * rbcExampleShockVals;
double * rbcExampleDataVals;
double * rbcExampleFP;
/*double * rbcExampleIntercept;*/
/*double * rbcExampleEasyPathQ;*/
/*double * rbcExampleTargetPathQ;*/
double * rbcExamplePathQ;
/*double * rbcExampleZeroPathQ;*/
