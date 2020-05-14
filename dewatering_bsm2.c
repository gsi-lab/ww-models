/*
 * dewatering_bsm2.c calculates the water and sludge stream concentrations   
 * from an 'ideal' dewatering unit based on a fixed percentage of solids in
 * the dewatered sludge. A defined amount of total solids are removed from
 * the influent sludge stream and goes into the stream of dewatered sludge
 * and the remaining will leave with the reject water phase. 
 * Soluble concentrations are not affected.
 * Temperature is also handled ideally, i.e. T(out)=T(in).
 *  
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 */

#define S_FUNCTION_NAME dewatering_bsm2

#include "simstruc.h"

#define PAR	ssGetArg(S,0)

int i;
/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 23);  /* number of inputs                      */
    ssSetNumOutputs(       S, 46);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 1);   /* number of input arguments             */
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
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double dewater_perc, X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS, TSS_removal_perc, dewater_factor, TSSin, Qu_factor, reject_factor;

  dewater_perc = mxGetPr(PAR)[0];
  TSS_removal_perc = mxGetPr(PAR)[1];
  X_I2TSS = mxGetPr(PAR)[2];
  X_S2TSS = mxGetPr(PAR)[3];
  X_BH2TSS = mxGetPr(PAR)[4];
  X_BA2TSS = mxGetPr(PAR)[5];
  X_P2TSS = mxGetPr(PAR)[6];
  
  TSSin = X_I2TSS*u[2]+X_S2TSS*u[3]+X_BH2TSS*u[4]+X_BA2TSS*u[5]+X_P2TSS*u[6]+X_BA2TSS*(u[20]+u[21]+u[22]);  /* cellulose is incorporated */
  dewater_factor = dewater_perc*10000.0/TSSin; 
  Qu_factor = TSS_removal_perc/(100.0*dewater_factor);
  reject_factor = (1.0-TSS_removal_perc/100.0)/(1.0-Qu_factor);
  
  if (dewater_factor > 1) {
    /* sludge */
    y[0]=u[0];
    y[1]=u[1];
    y[2]=u[2]*dewater_factor;
    y[3]=u[3]*dewater_factor;
    y[4]=u[4]*dewater_factor;
    y[5]=u[5]*dewater_factor;
    y[6]=u[6]*dewater_factor;
    y[7]=u[7];
    y[8]=u[8];
    y[9]=u[9];
    y[10]=u[10];
    y[11]=u[11]*dewater_factor;
    y[12]=u[12];
    y[13]=TSSin*dewater_factor;
    y[14]=u[14]*Qu_factor;

    y[15]=u[15]; /* Temp */
    
    /* Dummy states */
    y[16]=u[16];
    y[17]=u[17];
    y[18]=u[18];
    y[19]=u[19];
    y[20]=u[20]*dewater_factor;
    y[21]=u[21]*dewater_factor;
    y[22]=u[22]*dewater_factor; /* Cellulose */
   
    /* reject */
    y[23]=u[0];
    y[24]=u[1];
    y[25]=u[2]*reject_factor;
    y[26]=u[3]*reject_factor;
    y[27]=u[4]*reject_factor;
    y[28]=u[5]*reject_factor;
    y[29]=u[6]*reject_factor;
    y[30]=u[7];
    y[31]=u[8];
    y[32]=u[9];
    y[33]=u[10];
    y[34]=u[11]*reject_factor;
    y[35]=u[12];
    y[36]=TSSin*reject_factor;
    y[37]=u[14]*(1.0-Qu_factor);
    
    y[38]=u[15]; /* Temp */
    
    /* Dummy states */
    y[39]=u[16];
    y[40]=u[17];
    y[41]=u[18];
    y[42]=u[19];
    y[43]=u[20]*reject_factor;  
    y[44]=u[21]*reject_factor;
    y[45]=u[22]*reject_factor;
    }
  else 
  /* the influent is too high on solids to thicken further */
  /* all the influent leaves with the dewatered flow */
  {
    /* for sludge */
    for (i = 0; i < 23; i++) {
  
     y[i]= u[i];  
        
      }    
     y[13]=TSSin;
      
    /* for reject water */
     
     for (i = 23; i < 46; i++) {
  
     y[i]= 0;  
        
      }    
     
  }
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

