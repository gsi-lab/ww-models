/*
 * This model is developed for Enhance RBF based on the co-relation model used in Boichoi et al..,2019. 
 *  
 * Copyright: Chitta Ranjan Behera, CAPEC-PROCESS, DTU chemical, Denmark
 */

#define S_FUNCTION_NAME salsnes_bsm2N_withpoly_fastspeed

#include "simstruc.h"
#include <math.h> 

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
 	 ssSetNumContStates(    S, 0);     /* number of continuous states */
	 ssSetNumDiscStates(    S, 0);      /* number of discrete states */
	 ssSetNumInputs(        S, 23);     /* number of inputs  */
	 ssSetNumOutputs(       S, 46);     /* number of outputs */
	 ssSetDirectFeedThrough(S, 1);      /* direct feedthrough flag */
	 ssSetNumSampleTimes(   S, 1);      /* number of sample times */
	 ssSetNumInputArgs(     S, 0);      /* number of input arguments */
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

}


/*
 * mdlOutputs - compute the outputs
 */
static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
    double  Qu, Qo,f_ps,tss_eff,cod_par,cel_eff,XS_fra,XI_fra,XBH_fra,xnd_eff;
  
  double   X_S2TSS,X_BH2TSS,X_I2TSS;
  
  X_S2TSS = 1.01; // Boichoi et al..,2019
  X_BH2TSS = 0.35; // Boichoi et al..,2019
  X_I2TSS = 0.96; // Boichoi et al..,2019
  
  /* Calculation of efficiency and flow ratio Trojan model prediction */

  tss_eff=0.57;    /* Calculated separately running Trojan model with PI. 41 % for no poly, 49% for with poly 2 ppm, 57% fo with poly 4 ppm */
  XS_fra=tss_eff*X_S2TSS;
  XBH_fra=tss_eff*X_BH2TSS;
  XI_fra=tss_eff*X_I2TSS;
  xnd_eff=XS_fra;
  cel_eff=0.0; /* as it is not part of this influent */

  f_ps=0.003084957511230;  /* Calculated separately running Trojan model*/
  
         
 
  Qu = f_ps*u[14];         /* underflow from Salsnes to AD */
  Qo = u[14] - Qu;        /* overflow from Salsnes to AS */ 
   

   
 /* Calculation of ASM1 state outputs */
 
 /* effluent/overflow */
 
    y[0]=u[0];                /* S_I */
    y[1]=u[1];               /* S_S */
    y[2]=u[2]*(1.0-XI_fra);       /* X_I */
    y[3]=u[3]* (1.0-XS_fra);       /* X_S */
    y[4]=u[4]* (1.0-XBH_fra);    /* X_BH  */
    y[5]=u[5]* (1.0-XBH_fra);    /* X_BA  */
    y[6]=u[6]* (1.0-XI_fra);  /* X_P ,  */
    y[7]=u[7];                /* S_O */
    y[8]=u[8];                /* S_NO3 */
    y[9]=u[9];                /* S_NH */
    y[10]=u[10];              /* S_ND */
    y[11]=u[11]*(1.0-xnd_eff);  /* X_ND */ 
    y[12]=u[12];               /* S_ALK */
    y[13]=u[13]*(1.0-tss_eff); /* S_TSS, This TSS estimation is not correct one. */ 
    y[14]=Qo;
    y[15]=u[15];                 /* T */
    y[16]=u[16];                 /* S_NO2 */
    y[17]=u[17];                /* S_NO */
    y[18]=u[18];               /* S_N2O */
    y[19]=u[19];                   /* S_N2 */
    y[20]=u[20]* (1.0-XBH_fra);    /* X_BA2 */
    y[21]=u[21]* (1.0-XBH_fra);    /* X_ANAOB */                 
    y[22]=u[22]*(1.0-cel_eff); /* X_C */ 

   
        
 /* underflow/ sludge */
 
    y[23]=u[0];
    y[24]=u[1];
    y[25]=(u[14]*u[2]- Qo*y[2])/Qu; /* X_I */
    y[26]=(u[14]*u[3]- Qo*y[3])/Qu; /* X_S */
    y[27]=(u[14]*u[4]- Qo*y[4])/Qu; 
    y[28]=(u[14]*u[5]- Qo*y[5])/Qu; 
    y[29]=(u[14]*u[6]- Qo*y[6])/Qu; 
    y[30]=u[7];
    y[31]=u[8];
    y[32]=u[9];
    y[33]=u[10];
    y[34]=(u[14]*u[11]- Qo*y[11])/Qu; /* X_ND */
    y[35]=u[12];
    y[36]=(u[14]*u[13]- Qo*y[13])/Qu; /* This TSS estimation is not correct when cellulose is consider. */
    y[37]=Qu;
    y[38]=u[15]; /* Temp */

    y[39]=u[16];
    y[40]=u[17];
    y[41]=u[18];
    y[42]=u[19];
    y[43]=(u[14]*u[20]- Qo*y[20])/Qu; 
    y[44]=(u[14]*u[21]- Qo*y[21])/Qu; 
    y[45]=(u[14]*u[22]- Qo*y[22])/Qu; /* Cellulose */
 
    
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

