
/*This code is developed for IFAS process. for 10 layers
 * in this version the the calculation for Uf is modified and tested.
 *
 * Chitta Ranjan Behera
 * DTU,  2017
 *NOTE: always check whther parameters values are properly passed or not! Ex- Temperature, density etc. 
*/ 

#define S_FUNCTION_NAME ifas_model_v3_10layers

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)
#define REACTOR	ssGetArg(S,1)
#define KIN_PAR	ssGetArg(S,2)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 204); /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 23);   /*   number of inputs                    */
    ssSetNumOutputs(       S, 206);   /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 3);   /* number of input arguments             */
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

for (i = 0; i < 204; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x ,double *u, SimStruct *S, int tid)
{
double L, V_bulk,V_reactor, z_max;
double xtemp[200], xtemp1[20][10];
double Q_in, Q_out, SO2ref, Temp_op;
double H_N2_base,H_N2_0,T_op,P_N2_air,H_N2_1,kLa_N2,kLa,D_N2,D_O2,N2_gas,H_N2;

int i,j;
double fill_frac;

fill_frac=mxGetPr(REACTOR)[5];

V_reactor=              x[203];
Q_in=       u[20];      
Q_out=      Q_in;           
  
z_max =                 x[201];     

V_bulk = V_reactor - 500.0*fill_frac*V_reactor*z_max;  /*surface area is constant, 40% filling*/


for (i = 0; i < 200; i++) {
    if (x[i] < 1e-12)
      xtemp[i] = 0.0;
   else
      xtemp[i] = x[i];  
}

for (i = 0; i < 20; i++) {
    for (j = 0; j < 10; j++) {
        xtemp1[i][j] = xtemp[j+i*10];
    }
 }

for (i = 0; i < 20; i++) {
    for (j = 0; j < 9; j++) {
    y[j+i*10] = xtemp1[i][j]; 
    }
 }

for (i = 0; i < 20; i++) {
    y[9+i*10] = xtemp[9+i*10]/V_bulk; 
}    

y[200] = 0.0;
y[201] = x[201];
y[202] = V_bulk ;
y[203] = x[203];
y[204] = Q_out;

// Temp_op = 15.0;
Temp_op = u[21];   //unlock it for BSM2N
T_op = 273.0 + Temp_op;
kLa= u[22];
D_O2=                   2.2e-4;                  
D_N2=                   1.6e-4;       

/* mass transfer parameters (N2) */
H_N2_base =             6.5E-4;
H_N2_0 =                H_N2_base*exp(1300.0*(1.0/293.0 - 1.0/T_op));
// H_N2_0 =                H_N2_base*exp(1300.0*(1.0/T_op - 1.0/293.0));
P_N2_air =              0.78;
H_N2_1  =               H_N2_0*1000.0*28.0;              
kLa_N2 =                1.0*kLa* sqrt(D_N2) / sqrt(D_O2);

 N2_gas = -kLa_N2  *(H_N2_1*P_N2_air    - y[49])*V_bulk/1000.0;  // Kg/d


y[205] = N2_gas;


}

/*
 * mdlUpdate - perform action at major integration time step
 */

static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}

/*
 * mdlDerivatives - compute the derivatives  !!!!! 
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)

{

double V_reactor ,Q_in, kLa, Q_out, n_granules,  
       pKa_amm, pKa_nitri, pKa_HCO3;
double theta, density, z_max;
double D_NH4, D_NO2, D_O2, D_NO3, D_N2, D_S, f, L_B, z_init, k_NH4, k_NO2, k_O2, k_NO3, k_N2, k_S, S_O2_sat, D_S_I, k_S_I, D_IC, k_IC,D_IP, k_IP, D_ION, k_ION, D_CO2, k_CO2, H_CO2_1, P_CO2_air, kLa_CO2, D_CH4;
double mu_max_AOB, K_O2_AOB, K_NH3_AOB, b_AOB, Y_AOB, mu_max_NOB, K_O2_NOB, K_HNO2_NOB, b_NOB, Y_NOB, 
       mu_max_AnAOB, K_O2_AnAOB, K_NH3_AnAOB, K_HNO2_AnAOB, b_AnAOB, Y_AnAOB, f_i, i_NXB, i_NXI, k_H, K_X, 
       mu_max_HB, K_O2_HB, K_NH4_HB, K_NO2_HB, K_NO3_HB, K_S_HB, nu_HB, b_HB, Y_HB, K_HNO2_AOB,
       C_SI, C_SS, C_XI, C_XS, C_XB, C_XP, P_SI, P_SS, P_XI, P_XS, P_XB, P_XP, C_CH4;
double S_NH4_in, S_O2_in, S_NO2_in, S_NO3_in, S_N2_in, S_S_in, S_I_in, S_IC_in, S_IP_in, X_AOB_in, X_NOB_in, X_AnAOB_in, X_HB_in, X_S_in, X_I_in;

double diffus[10];
double transfer[20];
double influent[20];
double cc[204];

double c[20][10];
double xtemp[200];

double C_NH3[10],C_HNO2[10], C_HCO3[10],C_bulk[20], pH;
double C_NH3_bulk, C_HNO2_bulk, C_CO2_bulk, pH_bulk;

double r[20][10];
double r_bulk[20];

double ru_VSS[10], ru_ISS[10];

double sum;
double L, delta_x, A_biofilm, trapz_area;
double z[9];
double dudz[9];
double ru [9];
double A[9];
double product[9];
double uF[9];
double uD, uB, V_bulk;
double dCdt[20][10];
double dCdx2;
double dLdt;
double dVdt;
double ION_Na_in, ION_K_in, ION_Cl_in, ION_Ca_in, ION_Mg_in, ION_SO4_in, ION_Al_in, ION_Fe2_in, ION_Fe3_in, ION_HS_in, ION_NO2_in, S_CH4_in;
double stripping_CO2;
double K_ALK_AOB, K_pH_AOB, pHopt_AOB, K_pH_NOB, pHopt_NOB, K_pH_AnAOB, pHopt_AnAOB;
double N_SI, N_SS, N_XI, N_XS, N_XB, N_XP;
double H_N2_1, P_N2_air, kLa_N2;
double D_ionNa, D_ionK, D_ionCl, D_ionCa, D_ionMg, D_ionSO4, D_ionAl, D_ionFe2, D_ionFe3, D_ionHS; 

double T_base,T_op, R_ct, H_CO2_base,H_CO2_0, factor, H_N2_base, H_N2_0  ;

double H_H2S_base,H_H2S_0, P_H2S_air, H_H2S_1, kLa_H2S; 
double H_NH3_base,H_NH3_0, P_NH3_air, H_NH3_1, kLa_NH3; 
double H_CH4_base,H_CH4_0, P_CH4_air, H_CH4_1, kLa_CH4; 


// NEW

double stoich0_n, stoich1_n, stoich2_n, stoich3_n, stoich4_n, stoich5_n, stoich6_n, stoich7_n, stoich8_n, stoich9_n, stoich10_n, stoich11_n, stoich12_n, stoich13_n, stoich14_n, stoich15_n, stoich16_n, stoich17_n, stoich18_n, stoich19_n, stoich20_n, stoich21_n, stoich22_n, stoich23_n, stoich24_n, stoich25_n, stoich26_n, stoich27_n, stoich28_n, stoich29_n, stoich30_n, stoich31_n, stoich32_n, stoich33_n, stoich34_n, stoich35_n, stoich36_n, stoich37_n, stoich38_n, stoich39_n, stoich40_n; // NEW
double proc_bulk[15], reac_bulk[20];
double proc_bio[15][10], reac_bio[20][10];
double time;

double pH_set,Temp_op,fill_frac;

double eps;
double R_struP_bulk, R_CalcP_bulk,R_ACPP_bulk ;
double P1_I_in, P2_I_in, P3_I_in, P4_I_in;


double COD_TNN, COD_NO3, COD_NO2, COD_N2;
double TAN_TNN,TNN_NO3,TNN_N2,NO3_N2,NO3_NO2,NO2_N2;
double COD_SO4,COD_HS,COD_S0;
double HS_S0,S0_SO4,HS_SO4;
double NO3N2_HSS0,NO2N2_HSS0,NO3NO2_S0SO4;
double pKa_H2S, C_H2S_bulk, pKa_CO2;

/* i = current species
 * j = current discretized point in biofilm
 * k = current process
 */
int i, j, k, nx;

// Number of layers in biofilm

nx = 10.0;

// biofilm parameter
      
density =   mxGetPr(REACTOR)[1];    
z_max =     mxGetPr(REACTOR)[2];
pH_set =  mxGetPr(REACTOR)[3];  // pH
S_O2_sat=   mxGetPr(REACTOR)[4];
fill_frac=mxGetPr(REACTOR)[5];

Temp_op=   u[21];  //unlock it when running for BSM2N
// Temp_op = 15.0;
// Design and operational parameters

Q_in=                   u[20];   
Q_out= Q_in;
kLa=                    u[22]; //??? check

T_base=                 273.0 + 20.0;
T_op =                  273.0 + Temp_op;        
R_ct =                  0.083145;        
V_reactor=              x[203];

////////////////////////////////////////////////////////////////////

COD_TNN = (-24.0-16.0*2.0+8.0)/14.0;
COD_NO3 = (-24.0-16.0*3.0+8.0)/14.0;
COD_NO2 = (-24.0-16.0*2.0+8.0)/14.0;
COD_N2 =  -48.0/28.0;
TAN_TNN = abs(0.0-COD_TNN);
TNN_NO3 = abs(COD_TNN-COD_NO3);
TNN_N2 =  abs(COD_TNN-COD_N2);
NO3_N2 =  abs(COD_NO3-COD_N2);
NO3_NO2 = abs(COD_NO3-COD_NO2);
NO2_N2  = abs(COD_NO3-COD_N2);
COD_SO4 = 0.0;
COD_HS = 64.0/32.0;
COD_S0 = 48.0/32.0;
HS_S0 = abs(COD_HS-COD_S0);
S0_SO4 = abs(COD_S0-COD_SO4);
HS_SO4 = abs(COD_HS-COD_SO4);
NO3N2_HSS0 = NO3_N2/HS_S0;
NO2N2_HSS0 = NO2_N2/HS_S0;
NO3NO2_S0SO4 = NO3_NO2/S0_SO4;

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// kinetic and stochiometric parameters (biochemical)
mu_max_AOB = (mxGetPr(KIN_PAR)[0])*exp(-0.094*(T_base - T_op));                        // [1/d]
K_O2_AOB = mxGetPr(KIN_PAR)[1];                                                    // [gCOD/m3]
K_NH3_AOB = mxGetPr(KIN_PAR)[2];                                                // [gN/m3]
b_AOB = (mxGetPr(KIN_PAR)[3])*exp(-0.094*(T_base - T_op));                           // [1/d]
Y_AOB = mxGetPr(KIN_PAR)[4];                                                     // [gCOD/gN]
mu_max_NOB = (mxGetPr(KIN_PAR)[5])*exp(-0.061*(T_base - T_op));                      // [1/d]
K_O2_NOB = mxGetPr(KIN_PAR)[6];                                                   // [gCOD/m3]
K_HNO2_NOB = mxGetPr(KIN_PAR)[7];                                            // [gN/m3]
b_NOB = (mxGetPr(KIN_PAR)[8])*exp(-0.061*(T_base - T_op));                           // [1/d]
Y_NOB = mxGetPr(KIN_PAR)[9];                                                    // [gCOD/gN]
mu_max_AnAOB = (mxGetPr(KIN_PAR)[10])*exp(-0.096*(T_base - T_op));                   // [1/d]
K_O2_AnAOB = mxGetPr(KIN_PAR)[11];                                              // Inhibition coefficient[gCOD/m3]
K_NH3_AnAOB = mxGetPr(KIN_PAR)[12];                                            // [gN/m3]
K_HNO2_AnAOB = mxGetPr(KIN_PAR)[13];                                           // [gN/m3]
b_AnAOB = (mxGetPr(KIN_PAR)[14])*exp(-0.096*(T_base - T_op));                         // [1/d]
Y_AnAOB = mxGetPr(KIN_PAR)[15];                                                    // [gCOD/gN]
f_i = mxGetPr(KIN_PAR)[16];                                                   // fraction converted to inert biomass [gCOD/gCOD] 
i_NXB = mxGetPr(KIN_PAR)[17];
i_NXI = mxGetPr(KIN_PAR)[18];
k_H = (mxGetPr(KIN_PAR)[19])*exp(-0.110*(T_base - T_op));                             // hydrolysis rate [1/d]
K_X = (mxGetPr(KIN_PAR)[20])*exp(-0.110*(T_base - T_op));                             // hydrolysis half saturation constant [gCOD/gCOD]
mu_max_HB = (mxGetPr(KIN_PAR)[21])*exp(-0.069*(T_base - T_op));                        // [1/d]
K_O2_HB = mxGetPr(KIN_PAR)[22];                                                  // [gCOD/m3]
K_NH4_HB = mxGetPr(KIN_PAR)[23];                                                  // [gN/m3]
K_NO2_HB = mxGetPr(KIN_PAR)[24];                                                   // [gN/m3]
K_NO3_HB = mxGetPr(KIN_PAR)[25];                                                   // [gN/m3]
K_S_HB = mxGetPr(KIN_PAR)[26];                                                           // [gCOD/m3]
nu_HB = mxGetPr(KIN_PAR)[27];                                                  // anoxic correction factor
b_HB = (mxGetPr(KIN_PAR)[28])*exp(-0.069*(T_base - T_op));                            // [1/d]
Y_HB = mxGetPr(KIN_PAR)[29];                                                       // [gCOD/gN]
K_HNO2_AOB = mxGetPr(KIN_PAR)[30];                                                     // [gN/m3], inhibition constant, approximation 



K_pH_AOB=8.21;
pHopt_AOB=7.23;
K_pH_NOB=8.21;
pHopt_NOB=7.23;
K_pH_AnAOB=8.21;
pHopt_AnAOB=7.23;

pKa_amm=9.25;    
pKa_nitri=3.25;  
pKa_CO2 =6.35;  
pKa_H2S = 7.2;


//////////////////////////////////////////////////////////////////////////////

D_NH4=                  1.957E-9*3600.0*24.0;    
D_NO2=                  1.912E-9*3600.0*24.0;    
D_O2=                   2.2e-4;                 
D_NO3=                  1.902E-9*3600.0*24.0;   
D_N2=                   1.6e-4;                 
D_S=                    1.0e-4;                 
D_S_I=                  1.0e-4;                  
D_IC=                   2.586E-9*3600.0*24.0;         
D_IP=                   0.824E-9*3600.0*24.0;        

f=                      0.75;                   
L_B=                    1.0e-5;         



// C, N and P fractions of the biochemical state variables
C_SI =                  0.36178; 
C_SS =                  0.31843; 
C_XI =                  0.36178; 
C_XS =                  0.31843; 
C_XB =                  0.36612;
C_XP =                  0.36612;
C_CH4 =                 0.1875; 

N_SI    =               0.06003;     
N_SS    =               0.03352;     
N_XI    =               0.06003;     
N_XS    =               0.03352;    
N_XB    =               0.08615;   
N_XP    =               0.06003;


influent [0] =      u[0];
influent [1] =      u[1];
influent [2] =      u[2];
influent [3] =      u[3];
influent [4] =      u[4];
influent [5] =      u[5];
influent [6] =      u[6];
influent [7] =      u[7]; 
influent [8] =      u[8]; 
influent [9] =      u[9]; 
influent [10] =     u[10]; 
influent [11] =     u[11]; 
influent [12] =     u[12]; 
influent [13] =     u[13]; 
influent [14] =     u[14]; 
influent [15] =     u[15]; 
influent [16] =     u[16]; 
influent [17] =     u[17]; 
influent [18] =     u[18]; 
influent [19] =     u[19]; 

diffus [0] =        f*D_NH4;
diffus [1] =        f*D_O2;
diffus [2] =        f*D_NO2;
diffus [3] =        f*D_NO3;
diffus [4] =        f*D_N2;
diffus [5] =        f*D_S;
diffus [6] =        f*D_S_I;
diffus [7] =        f*D_IC;            
diffus [8] =        f*D_IP;            
      
diffus [9] =       0.0;


transfer [0] =      D_NH4/L_B;
transfer [1] =      D_O2/L_B;
transfer [2] =      D_NO2/L_B;
transfer [3] =      D_NO3/L_B;
transfer [4] =      D_N2/L_B;
transfer [5] =      D_S/L_B;
transfer [6] =      D_S_I/L_B;  
transfer [7] =      D_IC/L_B; 
transfer [8] =      D_IP/L_B; 
  
transfer [9] =     0.0;

// // mass transfer parameters (DO)
 
S_O2_sat=               S_O2_sat;


/////////////////////////////////////////////////////////////////////////////////
/// bulk and biofilm temporal variables definition

L = x[201];
delta_x = L/9.0;

V_bulk = V_reactor - 500.0*fill_frac*V_reactor*L;
A_biofilm = 500.0*fill_frac*V_reactor;/*total surface area available */


for (j = 0; j < 200; j++) {
   if (x[j] < 0.0)
     xtemp[j] =     0.0;
   else
     xtemp[j] =     x[j];
}

for (j = 0; j < 200; j++) {
     cc[j] =        xtemp[j];
}

for (i = 0; i < 20; i++) {
    for (j = 0; j < 10; j++) {
     c[i][j] =      cc[j+i*10];
    }
}

for (i = 0; i < 20; i++) {
     C_bulk[i] =   c[i][10-1]/V_bulk;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/// (bio)chemical kinetics calculation (bulk)  



pH_bulk=        pH_set;
C_CO2_bulk=     C_bulk[7] - C_bulk[7]/(1+pow(10,(-pH_bulk+pKa_CO2))); 
C_NH3_bulk =    C_bulk[0]/(1+pow(10,(-pH_bulk+pKa_amm)));
C_HNO2_bulk =   C_bulk[2]/(1+pow(10,(pH_bulk-pKa_nitri)));
C_H2S_bulk =    C_bulk[18]/(1+pow(10,(pH_bulk-pKa_H2S)));



proc_bulk[0] =  mu_max_AOB   *C_NH3_bulk/(C_NH3_bulk+K_NH3_AOB)   *K_HNO2_AOB/(C_HNO2_bulk+K_HNO2_AOB)   *C_bulk[1]/(C_bulk[1]+K_O2_AOB)      *C_bulk[10];
proc_bulk[1] =  mu_max_NOB   *C_HNO2_bulk/(C_HNO2_bulk+K_HNO2_NOB)*C_bulk[1]/(C_bulk[1]+K_O2_NOB)                                             *C_bulk[11];    
proc_bulk[2] =  mu_max_AnAOB *C_NH3_bulk/(C_NH3_bulk+K_NH3_AnAOB) *C_HNO2_bulk/(C_HNO2_bulk+K_HNO2_AnAOB)*K_O2_AnAOB/(C_bulk[1]+K_O2_AnAOB)   *C_bulk[12];    
proc_bulk[3] =  b_AOB  *C_bulk[10];
proc_bulk[4] =  b_NOB  *C_bulk[11];
proc_bulk[5] =  b_AnAOB*C_bulk[12];
proc_bulk[6] =  mu_max_HB                 *C_bulk[1]/(C_bulk[1]+K_O2_HB) *C_bulk[5]/(C_bulk[5]+K_S_HB)*C_bulk[0]/(C_bulk[0]+K_NH4_HB)                            *C_bulk[13];
proc_bulk[7] =  mu_max_HB*nu_HB           *C_bulk[2]/(C_bulk[2]+K_NO2_HB)*C_bulk[5]/(C_bulk[5]+K_S_HB)*C_bulk[0]/(C_bulk[0]+K_NH4_HB)*K_O2_HB/(C_bulk[1]+K_O2_HB)*C_bulk[13];
proc_bulk[8] =  mu_max_HB*nu_HB           *C_bulk[3]/(C_bulk[3]+K_NO3_HB)*C_bulk[5]/(C_bulk[5]+K_S_HB)*C_bulk[0]/(C_bulk[0]+K_NH4_HB)*K_O2_HB/(C_bulk[1]+K_O2_HB)*C_bulk[13];
proc_bulk[9] =  b_HB*C_bulk[13];
proc_bulk[10] = k_H*C_bulk[14]/K_X;
proc_bulk[11] = kLa     *(S_O2_sat           - C_bulk[1]);                  
proc_bulk[12] = 0.0;   
proc_bulk[13] = 0.0;
proc_bulk[14] = 0.0;
proc_bulk[15] = 0.0;

/*SNH */     reac_bulk[0] = (-1.0 / Y_AOB - N_XB)*proc_bulk[0] + (-N_XB)*proc_bulk[1] + (-1.0 / Y_AnAOB - N_XB)*proc_bulk[2] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bulk[3] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bulk[4] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bulk[5] + (-(N_XB - N_SS / Y_HB))*proc_bulk[6] + (-(N_XB - N_SS / Y_HB))*proc_bulk[7] + (-(N_XB - N_SS / Y_HB))*proc_bulk[8] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bulk[9] + (-(N_SS - N_XS))*proc_bulk[10];
/*SO2 */     reac_bulk[1] = -(3.4285714 - Y_AOB) / (Y_AOB)*proc_bulk[0] - (1.1428571 - Y_NOB) / (Y_NOB)*proc_bulk[1] - (1 - Y_HB) / Y_HB*proc_bulk[6] + proc_bulk[11];
/*SNO2*/     reac_bulk[2] = 1.0 / Y_AOB*proc_bulk[0] - 1.0 / Y_NOB*proc_bulk[1] + (-1.0 / Y_AnAOB - 1 / 1.14257)*proc_bulk[2] - (1.0 - Y_HB) / (1.71428*Y_HB)*proc_bulk[7];
/*SNO3*/     reac_bulk[3] = 1.0 / Y_NOB*proc_bulk[1] + 1.0 / 1.14257*proc_bulk[2] - (1.0 - Y_HB) / (2.85714*Y_HB)*proc_bulk[8];
/*SN2 */     reac_bulk[4] = 2.0 / Y_AnAOB*proc_bulk[2] + (1.0 - Y_HB) / (1.71428*Y_HB)*proc_bulk[7] + (1.0 - Y_HB) / (2.85714*Y_HB)*proc_bulk[8];   /* no Kla, when Do manually set  */
/*SS */      reac_bulk[5] = -1.0 / Y_HB*proc_bulk[6] - 1.0 / Y_HB*proc_bulk[7] - 1.0 / Y_HB*proc_bulk[8] + 1.0*proc_bulk[10];
/*SI */      reac_bulk[6] = 0.0;
/*SiC*/      reac_bulk[7] = 0.0;
/*SiP*/      reac_bulk[8] = 0.0;
/*SSS*/      reac_bulk[9] = 0.0;

/*XAOB*/     reac_bulk[10] = 1.0*proc_bulk[0] - 1.0*proc_bulk[3];
/*XNOB*/     reac_bulk[11] = 1.0*proc_bulk[1] - 1.0*proc_bulk[4];
/*XANM*/     reac_bulk[12] = 1.0*proc_bulk[2] - 1.0*proc_bulk[5];
/*XOHO*/     reac_bulk[13] = 1.0*proc_bulk[6] + 1.0*proc_bulk[7] + 1.0*proc_bulk[8] - 1.0*proc_bulk[9];
/*XS*/       reac_bulk[14] = (1.0 - f_i)*proc_bulk[3] + (1.0 - f_i)*proc_bulk[4] + (1.0 - f_i)*proc_bulk[5] + (1.0 - f_i)*proc_bulk[9] - proc_bulk[10];
/*XI*/       reac_bulk[15] = f_i*proc_bulk[3] + f_i*proc_bulk[4] + f_i*proc_bulk[5] + f_i*proc_bulk[9];
/*XXX*/      reac_bulk[16] = 0.0;
/*XXX*/      reac_bulk[17] = 0.0;
/*XXX*/      reac_bulk[18] = 0.0;
/*XXX*/      reac_bulk[19] = 0.0;


//////////////////////////////////////////////////////////////////////////////////////////////////////
/// (bio)chemical kinetics calculation (biofilm)  

for (j = 0; j < 10; j++) {


//pH[j]= pH_set; 
pH= pH_set; 
 
///
C_NH3[j] = (c[0][j])/(1+pow(10,(-pH+pKa_amm)));
C_HNO2[j] =(c[2][j])/(1+pow(10,(pH-pKa_nitri)));  


proc_bio[0][j] = mu_max_AOB*  C_NH3[j] / (C_NH3[j] + K_NH3_AOB)     *K_HNO2_AOB / (C_HNO2[j] + K_HNO2_AOB)* c[1][j] / (c[1][j] + K_O2_AOB)     *c[10][j];
proc_bio[1][j] = mu_max_NOB  *C_HNO2[j] / (C_HNO2[j] + K_HNO2_NOB)  *c[1][j] / (c[1][j] + K_O2_NOB)                                        *c[11][j];
proc_bio[2][j] = mu_max_AnAOB*C_NH3[j] / (C_NH3[j] + K_NH3_AnAOB)   *C_HNO2[j] / (C_HNO2[j] + K_HNO2_AnAOB)*K_O2_AnAOB / (c[1][j] + K_O2_AnAOB)*c[12][j];
proc_bio[3][j] = b_AOB*  c[10][j];
proc_bio[4][j] = b_NOB*  c[11][j];
proc_bio[5][j] = b_AnAOB*c[12][j];
proc_bio[6][j] = mu_max_HB      *c[1][j] / (c[1][j] + K_O2_HB) *c[5][j] / (c[5][j] + K_S_HB)*c[0][j] / (c[0][j] + K_NH4_HB)                          *c[13][j];
proc_bio[7][j] = mu_max_HB*nu_HB*c[2][j] / (c[2][j] + K_NO2_HB)*c[5][j] / (c[5][j] + K_S_HB)*c[0][j] / (c[0][j] + K_NH4_HB)*K_O2_HB / (c[1][j] + K_O2_HB)*c[13][j];
proc_bio[8][j] = mu_max_HB*nu_HB*c[3][j] / (c[3][j] + K_NO3_HB)*c[5][j] / (c[5][j] + K_S_HB)*c[0][j] / (c[0][j] + K_NH4_HB)*K_O2_HB / (c[1][j] + K_O2_HB)*c[13][j];
proc_bio[9][j] = b_HB*c[13][j];
proc_bio[10][j] = k_H*c[14][j] / K_X;
proc_bio[11][j] = 0.0;
proc_bio[12][j] = 0.0;
proc_bio[13][j] = 0.0;
proc_bio[14][j] = 0.0;
proc_bio[15][j] = 0.0;


/*SNH4 */ reac_bio[0][j] = (-1.0 / Y_AOB - N_XB)*proc_bio[0][j] + (-N_XB)*proc_bio[1][j] + (-1.0 / Y_AnAOB - N_XB)*proc_bio[2][j] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bio[3][j] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bio[4][j] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bio[5][j] + (-(N_XB - N_SS / Y_HB))*proc_bio[6][j] + (-(N_XB - N_SS / Y_HB))*proc_bio[7][j] + (-(N_XB - N_SS / Y_HB))*proc_bio[8][j] + (-((1.0 - f_i)*N_XS + f_i*N_XI - N_XB))*proc_bio[9][j] + (-(N_SS - N_XS))*proc_bio[10][j];
/*SO2 */  reac_bio[1][j] = -(3.4285714 - Y_AOB) / (Y_AOB)*proc_bio[0][j] - (1.1428571 - Y_NOB) / (Y_NOB)*proc_bio[1][j] - (1 - Y_HB) / Y_HB*proc_bio[6][j];
/*SNO2*/  reac_bio[2][j] = 1.0 / Y_AOB*proc_bio[0][j] - 1.0 / Y_NOB*proc_bio[1][j] + (-1.0 / Y_AnAOB - 1.0 / 1.14257)*proc_bio[2][j] - (1.0 - Y_HB) / (1.71428*Y_HB)*proc_bio[7][j];
/*SNO3*/  reac_bio[3][j] = 1.0 / Y_NOB*proc_bio[1][j] + 1.0 / 1.14257*proc_bio[2][j] - (1.0 - Y_HB) / (2.8571*Y_HB)*proc_bio[8][j];
/*SN2 */  reac_bio[4][j] = 2.0 / Y_AnAOB*proc_bio[2][j] + (1.0 - Y_HB) / (1.71428*Y_HB)*proc_bio[7][j] + (1.0 - Y_HB) / (2.8571*Y_HB)*proc_bio[8][j];
/*SS */   reac_bio[5][j] = -1.0 / Y_HB*proc_bio[6][j] - 1.0 / Y_HB*proc_bio[7][j] - 1.0 / Y_HB*proc_bio[8][j] + 1.0*proc_bio[10][j];
/*SI */   reac_bio[6][j] = 0.0;
/*SiC */  reac_bio[7][j] = 0.0;
/*SiP */  reac_bio[8][j] = 0.0;
/*SSS */  reac_bio[9][j] = 0.0;

/*XAOB */ reac_bio[10][j] = 1.0*proc_bio[0][j] - 1.0*proc_bio[3][j];
/*XNOB */ reac_bio[11][j] = 1.0*proc_bio[1][j] - 1.0*proc_bio[4][j];
/*XANM */ reac_bio[12][j] = 1.0*proc_bio[2][j] - 1.0*proc_bio[5][j];
/*XOHO */ reac_bio[13][j] = 1.0*proc_bio[6][j] + 1.0*proc_bio[7][j] + 1.0*proc_bio[8][j] - 1.0*proc_bio[9][j];
/*XS */   reac_bio[14][j] = (1.0 - f_i)*proc_bio[3][j] + (1.0 - f_i)*proc_bio[4][j] + (1.0 - f_i)*proc_bio[5][j] + (1.0 - f_i)*proc_bio[9][j] - proc_bio[10][j];
/*XI */   reac_bio[15][j] = f_i*proc_bio[3][j] + f_i*proc_bio[4][j] + f_i*proc_bio[5][j] + f_i*proc_bio[9][j];
/*XXX */  reac_bio[16][j] = 0.0;
/*XXX */  reac_bio[17][j] = 0.0;
/*XXX */  reac_bio[18][j] = 0.0;
/*XXX */  reac_bio[19][j] = 0.0;


 
}

//////////////////////////////////////////////////////////////////////////////////
// Biofilm growth
for (j = 0; j < 10-1; j++) {
        
    ru_VSS[j] = reac_bio[10][j]+reac_bio[11][j]+reac_bio[12][j]+reac_bio[13][j]+reac_bio[14][j]+reac_bio[15][j]+reac_bio[16][j]+reac_bio[17][j]+reac_bio[18][j]+reac_bio[19][j];  // sum of all organisms. Vector with nx-1 elements 
    ru_ISS[j] = 0.0;  // sum of all organisms. Vector with nx-1 elements
    
    ru[j] = ru_VSS[j] +  ru_ISS[j];
}

for (j = 0; j < 10-1; j++) {
    dudz[j] = ru[j]/density;
    uF[j]= delta_x*dudz[j];  /* As it is a flat sheet , area is equal at each discretized point. so no trapezoidal*/
}


uD = pow((L/z_max),2)*uF[nx-2];
uB = uF[nx-2]-uD;


////////////////////////////////////////////////////////////////////////////////////////////
/// mass balancing (solubles)

for (i = 0; i < 10; i++) {
   
	 /* At the carrier surface, using forward difference and full distance */

	 dCdt[i][0] = diffus[i] * (c[i][1] - c[i][0]) / pow(delta_x, 2) + reac_bio[i][0];

	 /* After the surface (i.e from second node ), using central difference and half distance */
	 for (j = 1; j < nx - 2; j++) {

		 dCdt[i][j] = diffus[i] * (c[i][j + 1] - 2.0*c[i][j] + c[i][j - 1]) / pow(delta_x, 2) + reac_bio[i][j];
	 }

	 /* calculation the concentration at the biofilm interface (SL) */
	 c[i][nx - 2] = ((diffus[i] * c[i][nx - 3] / delta_x) + transfer[i] * C_bulk[i]) / ((diffus[i] / delta_x) + transfer[i]);

	 /* at the interface, using backwad difference. Note: no rate expression (r[i][j]) as only mass transfer happening */
	 dCdt[i][nx - 2] = diffus[i] * (c[i][nx - 2] - 2.0*c[i][nx - 3] + c[i][nx - 4]) / pow(delta_x, 2);

	 /*bulk */
	 dCdt[i][nx - 1] = Q_in*influent[i] - Q_out*C_bulk[i] - transfer[i] * (C_bulk[i] - c[i][nx - 2])*A_biofilm + reac_bulk[i] * V_bulk;
}



/// mass balancing (particulates)   

for (i= 10; i < 20; i++) {
	dCdt[i][0] = -(c[i][0] * dudz[0]) - (uF[0] * (c[i][1] - c[i][0]) / delta_x) + reac_bio[i][0];
	 	for (j = 1; j < nx-2; j++) {
   
       dCdt[i][j] = -(c[i][j] * dudz[j]) - (uF[j] * (c[i][j] - c[i][j - 1]) / delta_x) + reac_bio[i][j];
		}
		dCdt[i][nx - 2] = -(c[i][nx - 2] * dudz[nx - 2]) - (uD*(c[i][nx - 2] - c[i][nx - 3]) / delta_x) + reac_bio[i][nx - 2];
		
		dCdt[i][nx - 1] = Q_in*influent[i] - Q_out*C_bulk[i] + uD*c[i][nx - 2] * A_biofilm + reac_bulk[i] * V_bulk;  
	}



// Mass balancing

for (i = 0; i < 20; i++) {
    for (j = 0; j < 10; j++) {
        dx[j+i*10] = dCdt[i][j];
    }
}

for (i = 16; i < 20; i++) {
    for (j = 0; j < 10; j++) {
        dx[j+i*10] = 0.0;
    }
}


for (i = 9; i < 10; i++) {
    for (j = 0; j < 10; j++) {
        dx[j+i*10] = 0.0;
    }
}


dLdt = uB;
dVdt = Q_in-Q_out;

dx[200]   = 0.0; 
dx[201]   = dLdt; 
dx[202]   = 0.0;
dx[203]   = dVdt; 
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
