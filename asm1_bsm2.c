/*
 * asm1_bsm2 is a C-file S-function for IAWQ AS Model No 1 with temperature 
 * dependencies of the kinetic parameters. In addition to the ASM1 states, TSS
 * and dummy states are included. TEMPMODEL defines how temperature changes
 * in the input affects the reactor temperature. Temperature dependency for 
 * oxygen saturation concentration and KLa has also been added in accordance
 * with BSM2 documentation.
 *
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt et al., 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden 
 */

#define S_FUNCTION_NAME asm1_bsm2

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)
#define PAR	ssGetArg(S,1)
#define V	ssGetArg(S,2)
#define SOSAT	ssGetArg(S,3)
#define TEMPMODEL  ssGetArg(S,4)
#define ACTIVATE  ssGetArg(S,5)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 27);  /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 28);  /* number of inputs                      */
    ssSetNumOutputs(       S, 30);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 6);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
int i;

for (i = 0; i < 27; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{

double proc1, proc2,proc2x1, proc2x2, proc2x3, proc2x4, proc3,proc3x1,proc3x2, proc4, proc5, proc5x1, proc5x2, proc6, proc7, proc8, proc9, proc10;
double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13, reac16, reac17, reac18, reac19, reac20, reac21, reac22, reac23, reac24, reac25;

double vol, SO_sat, SO_sat_temp, KLa_temp;

double a_KlaN2O,	b_A1,	b_A2,	b_H,	b_KlaN2O,	D_N2,	D_N2O,	D_NO,	D_O2,	F_BOD_COD,	f_P,	F_TSS_COD,	H_N2,	H_N2O,	H_NO,	H_O2,	i_X_B,	i_X_P,	KlaN2O_anoxic,	k_a	, K_FA ,	K_FNA, k_h,	K_I10FA,	K_I10FNA,	K_I3NO,	K_I4NO,	K_I5NO,	K_I9FA,	K_I9FNA,	K_N2O,	K_NO,	K_NO2,	K_NO3,	K_OA1,	K_OA2,	K_OH,	K_OH1,	K_OH2,	K_OH3,	K_OH4,	K_OH5,	K_S1,	K_S2,	K_S3,	K_S4,	K_S5,	K_X	,mu_A1,	mu_A2,	mu_H,	n_g2,	n_g3,	n_g4,	n_g5,	n_h,	n_Y,	pH,	P_N2O_air,	P_N2_air,	P_NO_air,	P_O2_air, Y_A1, Y_A2, Y_H;
double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
double b_Ratkowsky_mu_A1,	b_Ratkowsky_mu_A2,	b_Ratkowsky_mu_H,	c_Ratkowsky_mu_A1,	c_Ratkowsky_mu_A2,	c_Ratkowsky_mu_H,	Temp_Ref,	theta_b_A1,	theta_b_A2,	theta_b_H,	theta_kla,	theta_k_a,	theta_k_h,	Tmax_Ratkowsky_mu_A1,	Tmax_Ratkowsky_mu_A2,	Tmax_Ratkowsky_mu_H,	Tmin_Ratkowsky_mu_A1,	Tmin_Ratkowsky_mu_A2,	Tmin_Ratkowsky_mu_H;
double mu_H_Temp, mu_A1_Temp, mu_A2_Temp, b_H_Temp, b_A1_Temp, b_A2_Temp, k_h_Temp, k_a_Temp;
double K_SNH_aob1,	K_SNH_aob2,	K_SNO2_aob,	K_SNO_aob,	K_SO_aob1,	K_SO_aob2,	Y_aob,	Y_nob, proc11,proc12;

/* ADD NEW PARAMETERS FOR ANAOB*/
double Y_AnAOB, mu_max_AnAOB, mu_AnAOB_Temp, K_NH3_AnAOB, K_HNO2_AnAOB, K_O2_AnAOB, b_AnAOB, b_AnAOB_Temp; 

double S_FA, KB_2_KW, K_A, S_FNA;
double S_NOX;

double Kla_N2O, Kla_NO, Kla_N2;
double FluxN2O_gas, FluxN2_gas, FluxNO_gas;
double K_SO_AOBden1,  K_IO_AOBden1,  K_SO_AOBden2,  K_IO_AOBden2,  n_AOB,  n_Y_AOB, K_FNA_aob, K_FA_aob, Tin;

int i;
double tempmodel, activate;  
  
    
    
   

a_KlaN2O	=       mxGetPr(PAR)[0];
b_A1	=           mxGetPr(PAR)[1];
b_A2 =          	mxGetPr(PAR)[2];
b_H =           	mxGetPr(PAR)[3];
b_KlaN2O	=       mxGetPr(PAR)[4];
D_N2	=           mxGetPr(PAR)[5];
D_N2O	=           mxGetPr(PAR)[6];
D_NO	=           mxGetPr(PAR)[7];
D_O2	=           mxGetPr(PAR)[8];
F_BOD_COD	=       mxGetPr(PAR)[9];
f_P  =              mxGetPr(PAR)[10];
F_TSS_COD	=       mxGetPr(PAR)[11];
H_N2	=           mxGetPr(PAR)[12];
H_N2O	=           mxGetPr(PAR)[13];
H_NO	=           mxGetPr(PAR)[14];
H_O2	=           mxGetPr(PAR)[15];
i_X_B	=           mxGetPr(PAR)[16];
i_X_P	=           mxGetPr(PAR)[17];
KlaN2O_anoxic	=   mxGetPr(PAR)[18];
k_a               = mxGetPr(PAR)[19];
K_FA               =mxGetPr(PAR)[20];
K_FNA              =mxGetPr(PAR)[21];
k_h  =              mxGetPr(PAR)[22];
K_I10FA	 =          mxGetPr(PAR)[23];
K_I10FNA	=       mxGetPr(PAR)[24];
K_I3NO	=           mxGetPr(PAR)[25];
K_I4NO	=           mxGetPr(PAR)[26];
K_I5NO	=           mxGetPr(PAR)[27];
K_I9FA  =           mxGetPr(PAR)[28];
K_I9FNA	=           mxGetPr(PAR)[29];
K_N2O	=           mxGetPr(PAR)[30];
K_NO	=           mxGetPr(PAR)[31];
K_NO2	=           mxGetPr(PAR)[32];
K_NO3   =           mxGetPr(PAR)[33];
K_OA1	=           mxGetPr(PAR)[34];
K_OA2	=           mxGetPr(PAR)[35];
K_OH	=           mxGetPr(PAR)[36];
K_OH1	=           mxGetPr(PAR)[37];
K_OH2	=           mxGetPr(PAR)[38];
K_OH3	=           mxGetPr(PAR)[39];
K_OH4	=           mxGetPr(PAR)[40];
K_OH5	=           mxGetPr(PAR)[41];
K_S1	=           mxGetPr(PAR)[42];
K_S2	=           mxGetPr(PAR)[43];
K_S3	=           mxGetPr(PAR)[44];
K_S4	=           mxGetPr(PAR)[45];
K_S5	=           mxGetPr(PAR)[46];
K_X	=               mxGetPr(PAR)[47];
mu_A1	=           mxGetPr(PAR)[48];
mu_A2	=           mxGetPr(PAR)[49];
mu_H	=           mxGetPr(PAR)[50];
n_g2	=           mxGetPr(PAR)[51];
n_g3	=           mxGetPr(PAR)[52];
n_g4	=           mxGetPr(PAR)[53];
n_g5	=           mxGetPr(PAR)[54];
n_h	    =           mxGetPr(PAR)[55];
n_Y	    =           mxGetPr(PAR)[56];
pH	    =           mxGetPr(PAR)[57];
P_N2O_air =     	mxGetPr(PAR)[58];
P_N2_air  =         mxGetPr(PAR)[59];
P_NO_air  =       	mxGetPr(PAR)[60];
P_O2_air =      	mxGetPr(PAR)[61];
Y_A1	=           mxGetPr(PAR)[62];
Y_A2	=           mxGetPr(PAR)[63];
Y_H	=               mxGetPr(PAR)[64];

X_I2TSS =           mxGetPr(PAR)[65];
X_S2TSS =           mxGetPr(PAR)[66];
X_BH2TSS =          mxGetPr(PAR)[67];
X_BA2TSS =          mxGetPr(PAR)[68];
X_P2TSS =           mxGetPr(PAR)[69];

b_Ratkowsky_mu_A1=	mxGetPr(PAR)[70];
b_Ratkowsky_mu_A2=	mxGetPr(PAR)[71];
b_Ratkowsky_mu_H=	mxGetPr(PAR)[72];
c_Ratkowsky_mu_A1=	mxGetPr(PAR)[73];
c_Ratkowsky_mu_A2=	mxGetPr(PAR)[74];
c_Ratkowsky_mu_H=	mxGetPr(PAR)[75];
Temp_Ref	=       mxGetPr(PAR)[76];
theta_b_A1	=       mxGetPr(PAR)[77];
theta_b_A2	=       mxGetPr(PAR)[78];
theta_b_H	=       mxGetPr(PAR)[79];
theta_kla	=       mxGetPr(PAR)[80];
theta_k_a	=       mxGetPr(PAR)[81];
theta_k_h	=       mxGetPr(PAR)[82];

Tmax_Ratkowsky_mu_A1	=   mxGetPr(PAR)[83];
Tmax_Ratkowsky_mu_A2	=   mxGetPr(PAR)[84];
Tmax_Ratkowsky_mu_H	    =   mxGetPr(PAR)[85];
Tmin_Ratkowsky_mu_A1	=   mxGetPr(PAR)[86];
Tmin_Ratkowsky_mu_A2	=   mxGetPr(PAR)[87];
Tmin_Ratkowsky_mu_H	    =   mxGetPr(PAR)[88];

K_SNH_aob1 = 	mxGetPr(PAR)[89];
K_SNH_aob2= 	mxGetPr(PAR)[90];
K_SNO2_aob= 	mxGetPr(PAR)[91];
K_SNO_aob=      mxGetPr(PAR)[92];
K_SO_aob1=      mxGetPr(PAR)[93];
K_SO_aob2=      mxGetPr(PAR)[94];
Y_aob	=       mxGetPr(PAR)[95];
Y_nob	=       mxGetPr(PAR)[96];

K_SO_AOBden1 =  mxGetPr(PAR)[97];
K_IO_AOBden1 =  mxGetPr(PAR)[98];
K_SO_AOBden2 =  mxGetPr(PAR)[99];
K_IO_AOBden2 =  mxGetPr(PAR)[100];
n_AOB        =  mxGetPr(PAR)[101];
n_Y_AOB      =  mxGetPr(PAR)[102];
K_FNA_aob    =  mxGetPr(PAR)[103];
K_FA_aob     =  mxGetPr(PAR)[104];


Y_AnAOB      =  mxGetPr(PAR)[105];
mu_max_AnAOB =  mxGetPr(PAR)[106];
K_NH3_AnAOB  =  mxGetPr(PAR)[107];
K_HNO2_AnAOB =  mxGetPr(PAR)[108];
K_O2_AnAOB   =  mxGetPr(PAR)[109];
b_AnAOB      =  mxGetPr(PAR)[110];
        
pH	    =  mxGetPr(PAR)[57];
vol =   mxGetPr(V)[0];
Tin=u[15];



KB_2_KW = exp(6344.0/(273.15 + Tin));
S_FA = (x[9]*pow(10,pH))/(KB_2_KW + pow(10,pH));

   
/* calculation of free nitrous acid */
K_A = exp(-2300.0/(273.15+ Tin));
S_FNA = (x[16] * 1.0 / (1.0 + K_A * pow(10,pH)));

mu_H_Temp  =    pow((b_Ratkowsky_mu_H *  (u[15]- Tmin_Ratkowsky_mu_H) *  (1.0 - exp(c_Ratkowsky_mu_H *  (u[15]- Tmax_Ratkowsky_mu_H)))),2.0);
mu_A1_Temp =    pow((b_Ratkowsky_mu_A1 * (u[15]- Tmin_Ratkowsky_mu_A1) * (1.0 - exp(c_Ratkowsky_mu_A1 * (u[15]- Tmax_Ratkowsky_mu_A1)))),2.0);

 
  tempmodel = mxGetPr(TEMPMODEL)[0];
  activate = mxGetPr(ACTIVATE)[0];
  
  for (i = 0; i < 13; i++) {                    /* ASM state variables: SI, SS, XI, XS, XBH, XBA1, XU, SO, SNO3, SNH, SND, XND, SALK*/
      y[i] = x[i];
  }

  y[13] = 0.75*x[2]+0.75*x[3]+0.75*x[4]+0.75*x[5]+0.75*x[6]+0.75*(x[20]+x[21])+0.75*x[22];  /* XTSS*/
 
   
  y[14] = u[14];                                  /* Flow */

  if (tempmodel < 0.5)                            /* Temp */ 
     y[15] = Tin;                                  
  else 
     y[15] = Tin; 
         
  
  
  y[16] = x[16];                                   /* 2 step nitrification four step denitrification (additional) variables: SNO2, SNO, SN2O,SN2 and XBA2*/
  y[17] = x[17];    
  y[18] = x[18];    
  y[19] = x[19];    
  y[20] = x[20];    
 y[22] = x[22]; /*Celluose fraction as inert*/
  /* dummy states, only give outputs if ACTIVATE = 1 */
  if (activate > 0.5) {
      
  y[21] = x[21];
   proc2x3 =   mu_H_Temp * n_g4 * x[4]* (x[1] / (K_S4 + x[1])) *(x[17]/(K_NO +  x[17]  + (pow(x[17],2)/K_I4NO)))  * (K_OH4 /(K_OH4 + x[7]  ));
  proc2x4 =   mu_H_Temp * n_g5 * x[4]* (x[1] / (K_S5 + x[1])) *(x[18]/(K_N2O + x[18])) *(K_OH5 / (K_OH5 + x[7])) * (K_I5NO/(K_I5NO + x[17]));
  proc12 =    n_AOB * (mu_A1_Temp * x[5] * (x[17] / (K_SNO_aob +  x[17])) * (S_FA / (K_FA_aob + S_FA)) * (x[7]/(K_SO_AOBden2 + (1.0-2.0*pow((K_SO_AOBden2/K_IO_AOBden2),0.5)) * x[7] + (pow(x[7],2.0)/K_IO_AOBden2))));/*NO red to N2O*/
  proc2x2 =   mu_H_Temp * n_g3 * x[4]* (x[1] / (K_S3 + x[1])) *(x[16]/(K_NO2 + x[16])) *(K_OH3 / (K_OH3 + x[7])) * (K_I3NO/(K_I3NO + x[17]));
    
  y[23] = vol*((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3;/*production of N2O by HB*/
  y[24] = vol*( - (1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4; /*consumption of N2O by HB*/
  y[25] = vol*2.0/Y_A1*proc12;  /* N2O produce by AOB*/
  y[26] = vol*((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2;  /* anoxic HB growth on NO2-*/
 
      
  }
  else if (activate < 0.5) {
  
  y[21] = 0.0;
  
  proc2x2 =   mu_H_Temp * n_g3 * x[4]* (x[1] / (K_S3 + x[1])) *(x[16]/(K_NO2 + x[16])) *(K_OH3 / (K_OH3 + x[7])) * (K_I3NO/(K_I3NO + x[17]));
  proc2x3 =   mu_H_Temp * n_g4 * x[4]* (x[1] / (K_S4 + x[1])) *(x[17]/(K_NO +  x[17]  + (pow(x[17],2.0)/K_I4NO)))  * (K_OH4 /(K_OH4 + x[7]  ));
  proc2x4 =   mu_H_Temp * n_g5 * x[4]* (x[1] / (K_S5 + x[1])) *(x[18]/(K_N2O + x[18])) *(K_OH5 / (K_OH5 + x[7])) * (K_I5NO/(K_I5NO + x[17]));
  proc12 =    n_AOB * (mu_A1_Temp * x[5] * (x[17] / (K_SNO_aob +  x[17])) * (S_FA / (K_FA_aob + S_FA)) * (x[7]/(K_SO_AOBden2 + (1.0-2.0*pow((K_SO_AOBden2/K_IO_AOBden2),0.5)) * x[7] + (pow(x[7],2.0)/K_IO_AOBden2))));/*NO red to N2O*/

  y[23] = vol*((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3;/*production of N2O by HB*/
  y[24] = vol*( - (1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4; /*consumption of N2O by HB*/
  y[25] = vol*2.0/Y_A1*proc12;  /* N2O produce by AOB*/
  y[26] = vol*((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2;  /* anoxic HB growth on NO2-*/
 
  }
  
  KLa_temp = u[27]*pow(theta_kla, (Tin-Temp_Ref)); /* not x[15]?*/
  
  Kla_N2O = pow(D_N2O,0.5) / pow(D_O2,0.5) * KLa_temp;
  Kla_NO =  pow(D_NO,0.5) / pow(D_N2O,0.5) * Kla_N2O;
  Kla_N2 =  pow(D_N2,0.5) / pow(D_N2O,0.5) * Kla_N2O;
  
  FluxNO_gas =  - Kla_NO *  ((P_NO_air *  14.0 / H_NO) -  x[17]) * vol;
  FluxN2O_gas = - Kla_N2O * ((P_N2O_air * 28.0 / H_N2O) - x[18]) * vol;
  FluxN2_gas =  - Kla_N2 *  ((P_N2_air * 28.0 /  H_N2) -  x[19]) * vol;
  
/* calculation of free ammonia   */
KB_2_KW = exp(6344.0/(273.15 + Tin));
S_FA = (x[9]*pow(10,pH))/(KB_2_KW + pow(10,pH));

   
/* calculation of free nitrous acid */
K_A = exp(-2300.0/(273.15+ Tin));
S_FNA = (x[16] * 1.0 / (1.0 + K_A * pow(10,pH)));

  
  y[27] =FluxNO_gas;
  y[28] =FluxN2O_gas; //K_A; 
  y[29] =FluxN2_gas; //S_FNA ;
  
}

/*
 * mdlUpdate - perform action at major integration time step
 */

static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}

/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{


double proc1, proc2,proc2x1, proc2x2, proc2x3, proc2x4, proc3,proc3x1,proc3x2, proc4, proc5, proc5x1, proc5x2, proc6, proc7, proc8, proc9, proc10;
double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13, reac16, reac17, reac18, reac19, reac20, reac21, reac22, reac23, reac24, reac25;

double vol, SO_sat, SO_sat_temp, KLa_temp;

double a_KlaN2O,	b_A1,	b_A2,	b_H,	b_KlaN2O,	D_N2,	D_N2O,	D_NO,	D_O2,	F_BOD_COD,	f_P,	F_TSS_COD,	H_N2,	H_N2O,	H_NO,	H_O2,	i_X_B,	i_X_P,	KlaN2O_anoxic,	k_a	, K_FA ,	K_FNA, k_h,	K_I10FA,	K_I10FNA,	K_I3NO,	K_I4NO,	K_I5NO,	K_I9FA,	K_I9FNA,	K_N2O,	K_NO,	K_NO2,	K_NO3,	K_OA1,	K_OA2,	K_OH,	K_OH1,	K_OH2,	K_OH3,	K_OH4,	K_OH5,	K_S1,	K_S2,	K_S3,	K_S4,	K_S5,	K_X	,mu_A1,	mu_A2,	mu_H,	n_g2,	n_g3,	n_g4,	n_g5,	n_h,	n_Y,	pH,	P_N2O_air,	P_N2_air,	P_NO_air,	P_O2_air, Y_A1, Y_A2, Y_H;
double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
double b_Ratkowsky_mu_A1,	b_Ratkowsky_mu_A2,	b_Ratkowsky_mu_H,	c_Ratkowsky_mu_A1,	c_Ratkowsky_mu_A2,	c_Ratkowsky_mu_H,	Temp_Ref,	theta_b_A1,	theta_b_A2,	theta_b_H,	theta_kla,	theta_k_a,	theta_k_h,	Tmax_Ratkowsky_mu_A1,	Tmax_Ratkowsky_mu_A2,	Tmax_Ratkowsky_mu_H,	Tmin_Ratkowsky_mu_A1,	Tmin_Ratkowsky_mu_A2,	Tmin_Ratkowsky_mu_H;
double mu_H_Temp, mu_A1_Temp, mu_A2_Temp, b_H_Temp, b_A1_Temp, b_A2_Temp, k_h_Temp, k_a_Temp;
double K_SNH_aob1,	K_SNH_aob2,	K_SNO2_aob,	K_SNO_aob,	K_SO_aob1,	K_SO_aob2,	Y_aob,	Y_nob, proc11,proc12;

/* ADD NEW PARAMETERS FOR ANAOB*/
double Y_AnAOB, mu_max_AnAOB, mu_AnAOB_Temp, K_NH3_AnAOB, K_HNO2_AnAOB, K_O2_AnAOB, b_AnAOB, b_AnAOB_Temp; 

double S_FA, KB_2_KW, K_A, S_FNA;
double S_NOX;

double Kla_N2O, Kla_NO, Kla_N2;
double FluxN2O_gas, FluxN2_gas, FluxNO_gas;
double K_SO_AOBden1,  K_IO_AOBden1,  K_SO_AOBden2,  K_IO_AOBden2,  n_AOB,  n_Y_AOB, K_FNA_aob, K_FA_aob, Tin;
double ;


double xtemp[27];
double tempmodel;
int i;

a_KlaN2O	=       mxGetPr(PAR)[0];
b_A1	=           mxGetPr(PAR)[1];
b_A2 =          	mxGetPr(PAR)[2];
b_H =           	mxGetPr(PAR)[3];
b_KlaN2O	=       mxGetPr(PAR)[4];
D_N2	=           mxGetPr(PAR)[5];
D_N2O	=           mxGetPr(PAR)[6];
D_NO	=           mxGetPr(PAR)[7];
D_O2	=           mxGetPr(PAR)[8];
F_BOD_COD	=       mxGetPr(PAR)[9];
f_P  =              mxGetPr(PAR)[10];
F_TSS_COD	=       mxGetPr(PAR)[11];
H_N2	=           mxGetPr(PAR)[12];
H_N2O	=           mxGetPr(PAR)[13];
H_NO	=           mxGetPr(PAR)[14];
H_O2	=           mxGetPr(PAR)[15];
i_X_B	=           mxGetPr(PAR)[16];
i_X_P	=           mxGetPr(PAR)[17];
KlaN2O_anoxic	=   mxGetPr(PAR)[18];
k_a               = mxGetPr(PAR)[19];
K_FA               =mxGetPr(PAR)[20];
K_FNA              =mxGetPr(PAR)[21];
k_h  =              mxGetPr(PAR)[22];
K_I10FA	 =          mxGetPr(PAR)[23];
K_I10FNA	=       mxGetPr(PAR)[24];
K_I3NO	=           mxGetPr(PAR)[25];
K_I4NO	=           mxGetPr(PAR)[26];
K_I5NO	=           mxGetPr(PAR)[27];
K_I9FA  =           mxGetPr(PAR)[28];
K_I9FNA	=           mxGetPr(PAR)[29];
K_N2O	=           mxGetPr(PAR)[30];
K_NO	=           mxGetPr(PAR)[31];
K_NO2	=           mxGetPr(PAR)[32];
K_NO3   =           mxGetPr(PAR)[33];
K_OA1	=           mxGetPr(PAR)[34];
K_OA2	=           mxGetPr(PAR)[35];
K_OH	=           mxGetPr(PAR)[36];
K_OH1	=           mxGetPr(PAR)[37];
K_OH2	=           mxGetPr(PAR)[38];
K_OH3	=           mxGetPr(PAR)[39];
K_OH4	=           mxGetPr(PAR)[40];
K_OH5	=           mxGetPr(PAR)[41];
K_S1	=           mxGetPr(PAR)[42];
K_S2	=           mxGetPr(PAR)[43];
K_S3	=           mxGetPr(PAR)[44];
K_S4	=           mxGetPr(PAR)[45];
K_S5	=           mxGetPr(PAR)[46];
K_X	=               mxGetPr(PAR)[47];
mu_A1	=           mxGetPr(PAR)[48];
mu_A2	=           mxGetPr(PAR)[49];
mu_H	=           mxGetPr(PAR)[50];
n_g2	=           mxGetPr(PAR)[51];
n_g3	=           mxGetPr(PAR)[52];
n_g4	=           mxGetPr(PAR)[53];
n_g5	=           mxGetPr(PAR)[54];
n_h	    =           mxGetPr(PAR)[55];
n_Y	    =           mxGetPr(PAR)[56];
pH	    =           mxGetPr(PAR)[57];
P_N2O_air =     	mxGetPr(PAR)[58];
P_N2_air  =         mxGetPr(PAR)[59];
P_NO_air  =       	mxGetPr(PAR)[60];
P_O2_air =      	mxGetPr(PAR)[61];
Y_A1	=           mxGetPr(PAR)[62];
Y_A2	=           mxGetPr(PAR)[63];
Y_H	=               mxGetPr(PAR)[64];

X_I2TSS =           mxGetPr(PAR)[65];
X_S2TSS =           mxGetPr(PAR)[66];
X_BH2TSS =          mxGetPr(PAR)[67];
X_BA2TSS =          mxGetPr(PAR)[68];
X_P2TSS =           mxGetPr(PAR)[69];

b_Ratkowsky_mu_A1=	mxGetPr(PAR)[70];
b_Ratkowsky_mu_A2=	mxGetPr(PAR)[71];
b_Ratkowsky_mu_H=	mxGetPr(PAR)[72];
c_Ratkowsky_mu_A1=	mxGetPr(PAR)[73];
c_Ratkowsky_mu_A2=	mxGetPr(PAR)[74];
c_Ratkowsky_mu_H=	mxGetPr(PAR)[75];
Temp_Ref	=       mxGetPr(PAR)[76];
theta_b_A1	=       mxGetPr(PAR)[77];
theta_b_A2	=       mxGetPr(PAR)[78];
theta_b_H	=       mxGetPr(PAR)[79];
theta_kla	=       mxGetPr(PAR)[80];
theta_k_a	=       mxGetPr(PAR)[81];
theta_k_h	=       mxGetPr(PAR)[82];

Tmax_Ratkowsky_mu_A1	=   mxGetPr(PAR)[83];
Tmax_Ratkowsky_mu_A2	=   mxGetPr(PAR)[84];
Tmax_Ratkowsky_mu_H	    =   mxGetPr(PAR)[85];
Tmin_Ratkowsky_mu_A1	=   mxGetPr(PAR)[86];
Tmin_Ratkowsky_mu_A2	=   mxGetPr(PAR)[87];
Tmin_Ratkowsky_mu_H	    =   mxGetPr(PAR)[88];

K_SNH_aob1 = 	mxGetPr(PAR)[89];
K_SNH_aob2= 	mxGetPr(PAR)[90];
K_SNO2_aob= 	mxGetPr(PAR)[91];
K_SNO_aob=      mxGetPr(PAR)[92];
K_SO_aob1=      mxGetPr(PAR)[93];
K_SO_aob2=      mxGetPr(PAR)[94];
Y_aob	=       mxGetPr(PAR)[95];
Y_nob	=       mxGetPr(PAR)[96];

K_SO_AOBden1 =  mxGetPr(PAR)[97];
K_IO_AOBden1 =  mxGetPr(PAR)[98];
K_SO_AOBden2 =  mxGetPr(PAR)[99];
K_IO_AOBden2 =  mxGetPr(PAR)[100];
n_AOB        =  mxGetPr(PAR)[101];
n_Y_AOB      =  mxGetPr(PAR)[102];
K_FNA_aob    =  mxGetPr(PAR)[103];
K_FA_aob     =  mxGetPr(PAR)[104];


Y_AnAOB      =  mxGetPr(PAR)[105];
mu_max_AnAOB =  mxGetPr(PAR)[106];
K_NH3_AnAOB  =  mxGetPr(PAR)[107];
K_HNO2_AnAOB =  mxGetPr(PAR)[108];
K_O2_AnAOB   =  mxGetPr(PAR)[109];
b_AnAOB      =  mxGetPr(PAR)[110];
        
vol =               mxGetPr(V)[0];

SO_sat =            mxGetPr(SOSAT)[0];

tempmodel =         mxGetPr(TEMPMODEL)[0];
Tin=u[15];


/* temperature compensation */

mu_H_Temp  =    pow((b_Ratkowsky_mu_H *  (u[15]- Tmin_Ratkowsky_mu_H) *  (1.0 - exp(c_Ratkowsky_mu_H *  (u[15]- Tmax_Ratkowsky_mu_H)))),2.0);
mu_A1_Temp =    pow((b_Ratkowsky_mu_A1 * (u[15]- Tmin_Ratkowsky_mu_A1) * (1.0 - exp(c_Ratkowsky_mu_A1 * (u[15]- Tmax_Ratkowsky_mu_A1)))),2.0);
mu_A2_Temp =    pow((b_Ratkowsky_mu_A2 * (u[15]- Tmin_Ratkowsky_mu_A2) * (1.0 - exp(c_Ratkowsky_mu_A2 * (u[15]- Tmax_Ratkowsky_mu_A2)))),2.0);

mu_AnAOB_Temp = mu_max_AnAOB * exp(-0.096*(20.0-Tin));
b_H_Temp =      b_H  * pow(theta_b_H,Tin- Temp_Ref);
b_A1_Temp =     b_A1 * pow(theta_b_A1,Tin- Temp_Ref);
b_A2_Temp =     b_A2 * pow(theta_b_A2,Tin- Temp_Ref);
b_AnAOB_Temp  = b_AnAOB * exp(-0.096*(20.0- Tin)); /* TO BE CHECKED OUT*/
k_h_Temp =      k_h * pow(theta_k_h,Tin- Temp_Ref);
k_a_Temp =      k_a * pow(theta_k_a,Tin- Temp_Ref);

SO_sat_temp =   0.9997743214*   (8.0/  10.5*(  56.12*  6791.5*  exp(-  66.7354 + 87.4755/((  Tin+ 273.15)/  100.0) + 24.4526*  log((Tin+273.15)/   100.0)))); /* van't Hoff equation */
KLa_temp =      u[27]*pow(theta_kla, (Tin-Temp_Ref));


/* calculation of free ammonia   */
KB_2_KW = exp(6344.0/(273.15 + Tin));
S_FA = (x[9]*pow(10,pH))/(KB_2_KW + pow(10,pH));

   
/* calculation of free nitrous acid */
K_A = exp(-2300.0/(273.15+ Tin));
S_FNA = (x[16] * 1.0 / (1.0 + K_A * pow(10,pH)));


/* calculation of the total oxidazed nitrogen*/
S_NOX = x[8]+x[16]+x[17]+x[18];
   
/* specific KLa */
Kla_N2O = pow(D_N2O,0.5) / pow(D_O2,0.5) * KLa_temp;
Kla_NO =  pow(D_NO,0.5) / pow(D_N2O,0.5) * Kla_N2O;
Kla_N2 =  pow(D_N2,0.5) / pow(D_N2O,0.5) * Kla_N2O;

/* procesess */

for (i = 0; i < 27; i++) {
   if (x[i] < 0.0)
     xtemp[i] = 0.0;
   else
     xtemp[i] = x[i];
}

if (u[27] < 0.0)
      x[7] = fabs(u[27]);


proc1 =     mu_H_Temp *        xtemp[4]* (xtemp[1] / (K_S1 + xtemp[1])) *(xtemp[7]/ (K_OH1 + xtemp[7]));
proc2x1 =   mu_H_Temp * n_g2 * xtemp[4]* (xtemp[1] / (K_S2 + xtemp[1])) *(xtemp[8] /(K_NO3 + xtemp[8])) * (K_OH2 / (K_OH2 + xtemp[7]));
proc2x2 =   mu_H_Temp * n_g3 * xtemp[4]* (xtemp[1] / (K_S3 + xtemp[1])) *(xtemp[16]/(K_NO2 + xtemp[16])) *(K_OH3 / (K_OH3 + xtemp[7])) * (K_I3NO/(K_I3NO + xtemp[17]));
proc2x3 =   mu_H_Temp * n_g4 * xtemp[4]* (xtemp[1] / (K_S4 + xtemp[1])) *(xtemp[17]/(K_NO +  xtemp[17]  + (pow(xtemp[17],2.0)/K_I4NO)))  * (K_OH4 /(K_OH4 + xtemp[7]  ));
proc2x4 =   mu_H_Temp * n_g5 * xtemp[4]* (xtemp[1] / (K_S5 + xtemp[1])) *(xtemp[18]/(K_N2O + xtemp[18])) *(K_OH5 / (K_OH5 + xtemp[7])) * (K_I5NO/(K_I5NO + xtemp[17]));
proc3x1 =   mu_A1_Temp * xtemp[5]  * ( S_FA / (K_FA  + S_FA +  pow(S_FA,2)  / K_I9FA  )) * (xtemp[7] / (K_OA1 + xtemp[7])) * (K_I9FNA  / (K_I9FNA + S_FNA));
proc3x2 =   mu_A2_Temp * xtemp[20] * ( S_FNA /(K_FNA + S_FNA + pow(S_FNA,2) / K_I10FNA)) * (xtemp[7] / (K_OA2 + xtemp[7])) * (K_I10FA  / (K_I10FA + S_FA ));
proc4 =     b_H_Temp * xtemp[4];
proc5x1 =   b_A1_Temp * xtemp[5];
proc5x2 =   b_A2_Temp * xtemp[20];
proc6 =     k_a_Temp * xtemp[10]*xtemp[4];
proc7 =     k_h_Temp *(xtemp[3] /xtemp[4])/(K_X +(xtemp[3]/xtemp[4]))*((xtemp[7]/(K_OH+xtemp[7]))+n_h*(K_OH/(K_OH+xtemp[7]))*(S_NOX/(K_NO3+S_NOX)))*xtemp[4];
proc8 =     proc7*(xtemp[11]/xtemp[3]);
proc9 = mu_AnAOB_Temp*S_FA/(S_FA+K_NH3_AnAOB)*S_FNA/(S_FNA+K_HNO2_AnAOB)*K_O2_AnAOB/(xtemp[7]+K_O2_AnAOB)*xtemp[21];/*AnAOB growth*/
proc10 = b_AnAOB_Temp*xtemp[21];/*AnAOB decay*/
proc11 =     n_AOB * (mu_A1_Temp * xtemp[5] * (S_FNA / (K_FNA_aob + S_FNA)) * (S_FA / (K_FA_aob + S_FA)) * (xtemp[7]/(K_SO_AOBden1 + (1.0-2.0*pow((K_SO_AOBden1/K_IO_AOBden1),0.5)) * xtemp[7] + (pow(xtemp[7],2.0)/K_IO_AOBden1))));/*NO2 red to NO*/
proc12 =    n_AOB * (mu_A1_Temp * xtemp[5] * (xtemp[17] / (K_SNO_aob +  xtemp[17])) * (S_FA / (K_FA_aob + S_FA)) * (xtemp[7]/(K_SO_AOBden2 + (1.0-2.0*pow((K_SO_AOBden2/K_IO_AOBden2),0.5)) * xtemp[7] + (pow(xtemp[7],2.0)/K_IO_AOBden2))));/*NO red to N2O*/

/* reactions */

/* SI */    reac1 = 0.0;
/* SS */    reac2 = (- 1.0 / Y_H)*proc1 + ( - 1.0 / (Y_H *n_Y))*proc2x1 + (- 1.0 / (Y_H * n_Y))*proc2x2 + (- 1.0 / (Y_H * n_Y))*proc2x3 + (- 1.0 / (Y_H * n_Y))*proc2x4 + proc7;
/* XI */    reac3 = 0.0;
/* XS */    reac4 = (1.0 - f_P)*proc4 + (1.0 - f_P)*proc5x1 + (1.0 - f_P)*proc5x2 + (-1.0)*proc7 + proc10*(1.0-f_P);
/* XBH */   reac5 = (1.0)*proc1 + (1.0)*proc2x1 + (1.0)*proc2x2 + (1.0)*proc2x3 + (1.0)*proc2x4 + (-1.0)*proc4;
/* XBA1 */  reac6 =  proc3x1 - proc5x1 + proc11 + proc12;
/* XU */    reac7 = (f_P)*proc4 + (f_P)*proc5x1 + (f_P)*proc5x2+proc10*f_P;
/* SO2 */   reac8 =  - (1.0 - Y_H) / Y_H*proc1  - (3.4285714 - Y_A1) /Y_A1 *proc3x1  - (1.1428571 - Y_A2) / Y_A2*proc3x2 + (proc11+proc12)*(-(2.29-Y_A1)/Y_A1) ;
/* SNO3 */  reac9 = ( - (1.0 -Y_H * n_Y) / (1.1428571 * Y_H * n_Y))*proc2x1 + (1.0 / Y_A2)*proc3x2+proc9*1.52;
/* SNH  */  reac10 = ( - i_X_B)*proc1 + ( - i_X_B)*proc2x1 + ( - i_X_B)*proc2x2 + ( - i_X_B)*proc2x3 + ( - i_X_B)*proc2x4 + ( - i_X_B - (1.0 / Y_A1))*proc3x1 + ( - i_X_B)*proc3x2 + (1.0) * proc6 +proc9*(-1.0/Y_AnAOB-i_X_B) + (proc11+proc12)*(-i_X_B - 1.0/Y_A1);
/* SND  */  reac11 = -proc6+proc8;
/* XND */   reac12 = (i_X_B - f_P * i_X_P)*proc4 + (i_X_B - f_P * i_X_P)*proc5x1 + (i_X_B - f_P * i_X_P)*proc5x2 +(-1) *proc8+(i_X_B - f_P * i_X_P)*proc10;
/* SALK */  reac13 = ( - i_X_B / 14.0)*proc1 + (- i_X_B / 14.0)*proc2x1 + ( - (i_X_B / 14.0) + (1.0 - Y_H * n_Y) / (14.0 * ((3.4285714 - 2.8571429) * Y_H * n_Y)))*proc2x2 + ( - i_X_B / 14.0)*proc2x3 + ( - i_X_B / 14)*proc2x4 + (( - i_X_B) / 14.0 - 1.0 / (7.0 * Y_A1))*proc3x1 + ( - i_X_B / 14.0)*proc3x2 + (1.0 / 14.0)*proc6 + proc9*((-1/Y_AnAOB-i_X_B)/14.0-(-1.0/Y_AnAOB-1.52)/14.0-1.52/14.0)+proc11*((-i_X_B - 1.0/Y_A1)/14.0-(-1.0/Y_A1)/14.0)+proc12*((-i_X_B - 1.0/Y_A1)/14.0-(+1.0/Y_A1)/14.0);


/* SNO2 */reac16 =  ((1.0 - Y_H * n_Y) / (1.1428571 * Y_H * n_Y))*proc2x1 + ( - (1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2 + (1.0 / Y_A1)*proc3x1 + ( - 1.0 / Y_A2)*proc3x2  + proc9*(-1.0/Y_AnAOB-1.52)+proc11*(-1.0/Y_A1)+proc12*(1.0/Y_A1);
/* SNO  */reac17 =  ((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2 + ( - (1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3 + 2.0/Y_A1*proc11-2.0/Y_A1*proc12;
/* SN2O */reac18 =  ((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3 + ( - (1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4 + 2.0/Y_A1*proc12;
/* SN2  */reac19 =  ((1.0 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4+proc9*(2.0/Y_AnAOB);
/* XBA2 */reac20 = proc3x2 -proc5x2;

/*XAnAOB*/ reac21 = proc9-proc10;
reac22 = 0.0;
reac23 = 0.0;
reac24 = 0.0;
reac25 = 0.0;




/* mass balances */

dx[0] = 1.0/vol*(u[14]*(u[0]-x[0]))+reac1;                                                  /* SI */                                      
dx[1] = 1.0/vol*(u[14]*(u[1]-x[1]))+reac2;                                                  /* SS */     
dx[2] = 1.0/vol*(u[14]*(u[2]-x[2]))+reac3;                                                  /* XI */     
dx[3] = 1.0/vol*(u[14]*(u[3]-x[3]))+reac4;                                                  /* XS */    
dx[4] = 1.0/vol*(u[14]*(u[4]-x[4]))+reac5;                                                  /* XBH */    
dx[5] = 1.0/vol*(u[14]*(u[5]-x[5]))+reac6;                                                  /* XBA1 */    
dx[6] = 1.0/vol*(u[14]*(u[6]-x[6]))+reac7;                                                  /* XU */    
if (u[21] < 0.0)
      dx[7] = 0.0;
else
      dx[7] = 1.0/vol*(u[14]*(u[7]-x[7]))+reac8+KLa_temp*(SO_sat_temp - x[7]);

dx[8] = 1.0/vol*(u[14]*(u[8]-x[8]))+reac9;                                                  /* SNO3 */    
dx[9] = 1.0/vol*(u[14]*(u[9]-x[9]))+reac10;                                                 /* SNH */    
dx[10] = 1.0/vol*(u[14]*(u[10]-x[10]))+reac11;                                              /* SND */    
dx[11] = 1.0/vol*(u[14]*(u[11]-x[11]))+reac12;                                              /* XND */    
dx[12] = 1.0/vol*(u[14]*(u[12]-x[12]))+reac13;                                              /* SALK */    

dx[13] = 0.0;                                                                               /* TSS */

dx[14] = 0.0;                                                                               /* Flow */

if (tempmodel < 0.5)                                                                        /* Temp */    
   dx[15] = 0.0;                                  
else 
   dx[15] = 1.0/vol*(u[14]*(Tin-x[15]));  
//          dx[15] = 0.0;  
  
dx[16] = 1.0/vol*(u[14]*(u[16]-x[16]))+reac16;                                              /*SNO2 */

dx[17] = 1.0/vol*(u[14]*(u[17]-x[17]))+reac17 + Kla_NO * ((P_NO_air * 14.0 / H_NO) -  x[17]); /*SNO  */
dx[18] = 1.0/vol*(u[14]*(u[18]-x[18]))+reac18 + Kla_N2O *((P_N2O_air * 28.0 / H_N2O) - x[18]);/*SN2O */
dx[19] = 1.0/vol*(u[14]*(u[19]-x[19]))+reac19 + Kla_N2 * ((P_N2_air * 28.0 / H_N2) - x[19]);  /*SN2   */

dx[20] = 1.0/vol*(u[14]*(u[20]-x[20]))+reac20;                                              /*XBA2 */

dx[21] = 1.0/vol*(u[14]*(u[21]-x[21]))+reac21;                                              /*AnAOB */
dx[22] = 1.0/vol*(u[14]*(u[22]-x[22]));                                                     /*Cellulose */

dx[23] = 0.0;
dx[24] = 0.0;
dx[25] = 0.0;
dx[26] = 0.0;

//dx[23] = 1.0/vol*(u[14]*(u[23]-x[23]))+reac22;                                              /*D2 */
//dx[24] = 1.0/vol*(u[14]*(u[24]-x[24]))+reac23;                                              /*D3 */
//dx[25] = 1.0/vol*(u[14]*(u[25]-x[25]))+reac24;                                              /*D4 */
//dx[26] = 1.0/vol*(u[14]*(u[26]-x[26]))+reac25;                                              /*D5 */
}


/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


