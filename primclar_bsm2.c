/*
 * This function, primclar_bsm2.c, is a C-file S-function implementation of the Otterpohl/Freund
 * primary clarifier model. The implementation is to a large extent based on an
 * implementation of the Otterpohl/Freund model by Dr Jens Alex, IFAK, Magdeburg. 
 * In addition to ASM1 states, the clarifier will also pass on TSS, Q, temp and 5 dummy
 * states to effluent and underflow.
 * If TEMPMODEL=0 then T(out)=T(in), if TEMPMODEL=1 then T(out) is a first-order
 * equation based on the heat content of the influent, the reactor and outflow.
 *  
 *Modified by Chitta to include cellulose as state variables
 *
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 */

#define S_FUNCTION_NAME primclar_bsm2

#include "simstruc.h"
#include <math.h>

#define XINIT    ssGetSFcnParam(S,0)	/* initial state */
#define XVEKTOR  ssGetSFcnParam(S,1)	/* setteability of ASM1 components */
#define PAR_P ssGetSFcnParam(S,2)       /* special parameters for primary clarifier */
#define VOL	 ssGetSFcnParam(S,3)		/* Volume of primary clarifier */
#define PAR   ssGetSFcnParam(S,4)       /* parameters, needed for calculation of TSSeff and TSSu from ASM1 states */
#define TEMPMODEL  ssGetArg(S,5)        /* select the temperature model */

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
 	 ssSetNumContStates(    S, 23);     /* number of continuous states */
	 ssSetNumDiscStates(    S, 0);      /* number of discrete states */
	 ssSetNumInputs(        S, 23);     /* number of inputs */
	 ssSetNumOutputs(       S, 69);     /* number of outputs */
	 ssSetDirectFeedThrough(S, 1);      /* direct feedthrough flag */
	 ssSetNumSampleTimes(   S, 1);      /* number of sample times */
	 ssSetNumInputArgs(     S, 6);      /* number of input arguments */
	 ssSetNumRWork(         S, 0);      /* number of real work vector elements */
	 ssSetNumIWork(         S, 0);      /* number of integer work vector elements */
	 ssSetNumPWork(         S, 0);      /* number of pointer work vector elements */

}


/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S){
	 ssSetSampleTime(S, 0, 0.0);     /*  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME) */
	 ssSetOffsetTime(S, 0, 0.0);
}


/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
int i;

for (i = 0; i < 23; i++) {
    x0[i] = mxGetPr(XINIT)[i];
}
}


/*
 * mdlOutputs - compute the outputs
 */
static void mdlOutputs(double *y, const double *x, const double *u, SimStruct *S, int tid)
{
  double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
  int i;
  double tt, vol, rho, nCOD, nX, K, E, f_PS, Qu;
  double *pptr;
  double ff;
  double tempmodel;

  X_I2TSS = mxGetPr(PAR)[19];
  X_S2TSS = mxGetPr(PAR)[20];
  X_BH2TSS = mxGetPr(PAR)[21];
  X_BA2TSS = mxGetPr(PAR)[22];
  X_P2TSS = mxGetPr(PAR)[23];

  pptr = mxGetPr(XVEKTOR);   /* Parameters indicating settleability (0 or 1) */
  vol = mxGetPr(VOL)[0];     /* Primary clarifier volume */
  rho = mxGetPr(PAR_P)[0];   /* Settler efficiency correction */
  K = mxGetPr(PAR_P)[1];     /* Average CODpart/CODtot ratio */
  f_PS = mxGetPr(PAR_P)[3];  /* Ratio of primary sludge flow rate to the influent flow */
	
  tempmodel = mxGetPr(TEMPMODEL)[0];
  
  Qu = f_PS*u[14];         /* underflow from PC */
  E = u[14]/Qu;	           /* thickening factor, u[14] is the influent flow rate */
  tt = vol/(x[14]+0.001);  /* hydraulic retention time within primary clarifier */

  nCOD = rho*(2.88*K-0.118)*(1.45+6.15*log(tt*24.0*60.0));	/* Total COD removal efficiency in primary clarifier, in % */

  nX = nCOD/K;		/* nX = nCOD*nfak; Removal efficiency for particulate COD in %, since assumption that soluble COD is not removed*/

  if (nX > 100.0) { 
     nX = 100.0;
  }
  
  if (nX < 0.0) {
     nX = 0.0;
  }

  for (i = 0; i < 13; i++) {		    /* Calculation of ASM1 state outputs */
      ff = (1.0-pptr[i]*nX/100.0);
      y[i] = ff*x[i];	                /* effluent */
      if (y[i] < 0.0)
         y[i] = 0.0;
      y[i+23] = ((1.0-ff)*E + ff)*x[i];	/* primary sludge */
      if (y[i+23] < 0.0)
	     y[i+23] = 0.0;
  }
  
  for (i = 16; i < 23; i++) {		    /* Calculation of dummy state outputs */
      ff = (1.0-pptr[i]*nX/100.0);
      y[i] = ff*x[i];	                /* effluent */
      if (y[i] < 0.0)
         y[i] = 0.0;
      y[i+23] = ((1.0-ff)*E + ff)*x[i];	/* primary sludge */
      if (y[i+23] < 0.0)
	     y[i+23] = 0.0;
  }
  
  y[13] = 0.75*(y[2]+y[3]+y[4]+(y[5]+y[20]+y[21])+y[6]+y[22]);     /* TSS effluent,(Cellulose is incorporated) */ 
  y[14] = u[14]-Qu;                     /* flow rate effluent */
  
  if (tempmodel < 0.5)      /* Temperature effluent */                      
     y[15] = u[15];                                  
  else 
     y[15] = x[15];         
    
  y[36] = 0.75*(y[24]+y[25]+y[26]+(y[27]+y[42]+y[43])+y[28]+y[44]);     /* TSS primary sludge (Cellulose is incorporated)*/
  y[37] = Qu;              /* primary sludge flow rate */
  
  if (tempmodel < 0.5)      /* Temperature primary sludge */                      
     y[38] = u[15];                                  
  else 
     y[38] = x[15];         
  
  for (i = 0; i < 13; i++) {			/* ASM1 states */
      y[i+46] = x[i];	                        
      if (y[i+46] < 0.0)
	     y[i+46] = 0.0;
  }    
  
  y[59] = 0.75*(x[2]+x[3]+x[4]+(x[20]+x[5]+x[21])+x[6]+x[22]); /* TSS (Cellulose is incorporated) */
  
  y[60] = u[14];            /* Flow */
    
  if (tempmodel < 0.5)      /* Temp */                      
     y[61] = u[15];                                  
  else 
     y[61] = x[15];
  
  for (i = 16; i < 23; i++)      
      y[i+46] = x[i];
      if (y[i+46] < 0) {	/* Dummy states */		
         y[i+46] = 0;
  }
}


/*
 * mdlUpdate - perform action at major integration time step
 */
static void mdlUpdate(double *x, const double *u, SimStruct *S, int tid){

}


/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, const double *x, const double *u, SimStruct *S, int tid)
{
  int i;
  double t_m, vol;
  double tempmodel;
  
  t_m = mxGetPr(PAR_P)[2];
  
  vol = mxGetPr(VOL)[0];
  tempmodel = mxGetPr(TEMPMODEL)[0];
  
  for (i = 0; i < 13; i++) {                /* ASM1 states */
      dx[i] = 1.0/vol*(u[14]*(u[i]-x[i]));  /* mixing */
  }
  
  dx[13] = 0.0;                      /* TSS */
  
  dx[14] = (u[14]-x[14])/t_m;        /* Flow */
  
  if (tempmodel < 0.5)               /* Temp */    
     dx[15] = 0.0;                                  
  else 
     dx[15] = 1.0/vol*(u[14]*(u[15]-x[15]));  
  
  for (i = 16; i < 23; i++) {               /* Dummy states */
      dx[i] = 1.0/vol*(u[14]*(u[i]-x[i]));  /* mixing */
  }
    
}


/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{

}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"       /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"        /* Code generation registration function */
#endif



