/*
 * New version (no 3) of the ASM1 to ADM1 interface based on discussions
 * within the IWA TG BSM community during 2002-2006. Now also including charge
 * balancing and temperature dependency for applicable parameters.
 * Model parameters are defined in adm1init_bsm2.m
 * u is the input in ASM1 terminology + extra dummy states, 21 variables
 * plus one extra input = dynamic pH from the ADM1 system (needed for 
 * accurate charge balancing - also used the ADM1 to ASM1 interface).
 * If temperature control of AD is used then the operational temperature
 * of the ADM1 should also be an input rather than a defined parameter.
 * Temperature in the ADM1 and the ASM1 to ADM1 and the ADM1 to ASM1 
 * interfaces should be identical at every time instant.
 * Input vector:
 * u[0] : Si = soluble inert organic material (g COD/m3)
 * u[1] : Ss = readily biodegradable substrate (g COD/m3)
 * u[2] : Xi = particulate inert organic material (g COD/m3)
 * u[3] : Xs = slowly biodegradable substrate (g COD/m3)
 * u[4] : Xbh = active heterotrophic biomass (g COD/m3)
 * u[5] : Xba = active autotrophic biomass (g COD/m3)
 * u[6] : Xp = particulate product arising from biomass decay (g COD/m3)
 * u[7] : So = oxygen (g -COD/m3)
 * u[8] : Sno = nitrate and nitrite nitrogen (g N/m3)
 * u[9] : Snh = ammonia and ammonium nitrogen (g N/m3)
 * u[10] : Snd = soluble biogradable organic nitrogen (g N/m3)
 * u[11] : Xnd = particulate biogradable organic nitrogen (g N/m3)
 * u[12] : Salk = alkalinity (mole HCO3-/m3)
 * u[13] : TSS = total suspended solids (internal use) (mg SS/l)
 * u[14] : flow rate (m3/d)
 * u[15] : temperature (deg C)
 * u[16:20] : dummy states for future use
 * u[21] : pH in the anaerobic digester
 *
 * y is the output in ADM1 terminology + extra dummy states, 33 variables
 * y[0] : Ssu = monosacharides (kg COD/m3)
 * y[1] : Saa = amino acids (kg COD/m3)
 * y[2] : Sfa = long chain fatty acids (LCFA) (kg COD/m3)
 * y[3] : Sva = total valerate (kg COD/m3)
 * y[4] : Sbu = total butyrate (kg COD/m3)
 * y[5] : Spro = total propionate (kg COD/m3)
 * y[6] : Sac = total acetate (kg COD/m3)
 * y[7] : Sh2 = hydrogen gas (kg COD/m3)
 * y[8] : Sch4 = methane gas (kg COD/m3)
 * y[9] : Sic = inorganic carbon (kmole C/m3)
 * y[10] : Sin = inorganic nitrogen (kmole N/m3)
 * y[11] : Si = soluble inerts (kg COD/m3)
 * y[12] : Xc = composites (kg COD/m3)
 * y[13] : Xch = carbohydrates (kg COD/m3)
 * y[14] : Xpr = proteins (kg COD/m3)
 * y[15] : Xli = lipids (kg COD/m3)
 * y[16] : Xsu = sugar degraders (kg COD/m3)
 * y[17] : Xaa = amino acid degraders (kg COD/m3)
 * y[18] : Xfa = LCFA degraders (kg COD/m3)
 * y[19] : Xc4 = valerate and butyrate degraders (kg COD/m3)
 * y[20] : Xpro = propionate degraders (kg COD/m3)
 * y[21] : Xac = acetate degraders (kg COD/m3)
 * y[22] : Xh2 = hydrogen degraders (kg COD/m3)
 * y[23] : Xi = particulate inerts (kg COD/m3)
 * y[24] : scat+ = cations (metallic ions, strong base) (kmole/m3)
 * y[25] : san- = anions (metallic ions, strong acid) (kmole/m3)
 * y[26] : flow rate (m3/d)
 * y[27] : temperature (deg C)
 * y[28:32] : dummy states for future use
 *
 * ASM1 --> ADM1 conversion, version 3 for BSM2
 * Copyright: John Copp, Primodal Inc., Canada; Ulf Jeppsson, Lund
 *            University, Sweden; Damien Batstone, Univ of Queensland,
 *            Australia, Ingmar Nopens, Univ of Ghent, Belgium,
 *            Marie-Noelle Pons, Nancy, France, Peter Vanrolleghem,
 *            Univ. Laval, Canada, Jens Alex, IFAK, Germany and 
 *            Eveline Volcke, Univ of Ghent, Belgium.
 */

#define S_FUNCTION_NAME asm2adm_v3_bsm2_byRiccardo


#include "simstruc.h"
#include <math.h>

#define PAR  	ssGetArg(S,0)
#define min(a,b) (((a)<(b))?(a):(b))
/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 23);  /* number of inputs                      */
    ssSetNumOutputs(       S, 136);  /* number of outputs                     */
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
	double 	CODequiv, fnaa, fnxc, fnbac, fxni, fsni, fsni_adm, frlixs, frlibac, frxs_adm, fdegrade_adm, frxs_as, fdegrade_as;
    double  R, T_base, T_op, pK_w_base, pK_a_va_base, pK_a_bu_base, pK_a_pro_base, pK_a_ac_base, pK_a_co2_base, pK_a_IN_base;
    double  pH_adm, pK_w, pK_a_co2, pK_a_IN, alfa_va, alfa_bu, alfa_pro, alfa_ac, alfa_co2, alfa_IN, alfa_NH, alfa_alk, alfa_NO, factor;
    double	CODdemand, remaina, remainb, remainc, remaind, ScatminusSan;
	double 	sorgn, xorgn, xprtemp, xprtemp2, xlitemp, xlitemp2, xlitemp3, xchtemp, xchtemp2, xchtemp3;
    double CODequiv_NO3, CODequiv_NO2, CODequiv_NO, CODequiv_N2O,iN_B;
    double  biomass, biomass_nobio, biomass_bioN, remainCOD, inertX, xc, noninertX, inertS, utemp[23], utemp2[23];
    double  Y_H , deltaBIOMASS1 , deltaBIOMASS2, deltaBIOMASS3 , deltaBIOMASS3bis , CODdemand_remaining1, CODdemand_remaining2, CODdemand_remaining3, CODdemand_remaining4;
    double deltaNH1, deltaNH2, deltaNH3, deltaNH4,deltaNHtot;
    double NtotIN, Nremaining,O2,O2_consumedA,O2_consumedB,O2remainingA, O2remainingB, NHremainingA,NHremainingB,NHremainingC ,NHremainingD, NHremainingE, NHremainingF, NHremainingG, NHremainingH, NHremainingI, NHremainingL, XBHremainingA,XBHremainingB, XBHremainingC, XBHremainingD, XBHremainingE, XBHremainingF, XBHremainingG, XBHremainingH, XBHremainingI,XBHremainingL ;
    double n_Y,NO3,NO2,NO,N2O,NO3_consumedA,NO3_consumedB,NO3remainingA,NO3remainingB,NO2_consumedA,NO2_consumedB,NO2remainingA,NO2remainingB,NO_consumedA,N2O_consumedA,NO_consumedB,N2O_consumedB,NOremainingA,N2OremainingA,NOremainingB,N2OremainingB;
    double f_P,iN_XP;
	int	i;
    
    /* parameters defined in adm1init_bsm2.m, INTERFACEPAR */
    CODequiv = mxGetPr(PAR)[0];
    fnaa = mxGetPr(PAR)[1];
    fnxc = mxGetPr(PAR)[2];
    fnbac = mxGetPr(PAR)[3];
    fxni = mxGetPr(PAR)[4];
    fsni = mxGetPr(PAR)[5];
    fsni_adm = mxGetPr(PAR)[6];
    frlixs = mxGetPr(PAR)[7];
    frlibac = mxGetPr(PAR)[8];
    frxs_adm = mxGetPr(PAR)[9];
    fdegrade_adm = mxGetPr(PAR)[10];
    frxs_as = mxGetPr(PAR)[11];       /* not used in ASM2ADM */
    fdegrade_as = mxGetPr(PAR)[12];   /* not used in ASM2ADM */
    R = mxGetPr(PAR)[13]; 
    T_base = mxGetPr(PAR)[14];
    T_op = mxGetPr(PAR)[15];          /* should be an input variable if dynamic temperature control is used */
    pK_w_base = mxGetPr(PAR)[16];
    pK_a_va_base = mxGetPr(PAR)[17];
    pK_a_bu_base = mxGetPr(PAR)[18];
    pK_a_pro_base = mxGetPr(PAR)[19];
    pK_a_ac_base = mxGetPr(PAR)[20];
    pK_a_co2_base = mxGetPr(PAR)[21];
    pK_a_IN_base = mxGetPr(PAR)[22];  
    CODequiv_NO3 = mxGetPr(PAR)[23];
    CODequiv_NO2 = mxGetPr(PAR)[24];
    CODequiv_NO = mxGetPr(PAR)[25];
    CODequiv_N2O = mxGetPr(PAR)[26];
    iN_B=mxGetPr(PAR)[27];
    pH_adm = u[22];

    factor = (1.0/T_base - 1.0/T_op)/(100.0*R);
    pK_w = pK_w_base - log10(exp(55900.0*factor));
    pK_a_co2 = pK_a_co2_base - log10(exp(7646.0*factor));
    pK_a_IN = pK_a_IN_base - log10(exp(51965.0*factor));
    alfa_va = 1.0/208.0*(-1.0/(1.0 + pow(10, pK_a_va_base - pH_adm)));
    alfa_bu = 1.0/160.0*(-1.0/(1.0 + pow(10, pK_a_bu_base - pH_adm)));
    alfa_pro = 1.0/112.0*(-1.0/(1.0 + pow(10, pK_a_pro_base - pH_adm)));
    alfa_ac = 1.0/64.0*(-1.0/(1.0 + pow(10, pK_a_ac_base - pH_adm)));
    alfa_co2 = -1.0/(1.0 + pow(10, pK_a_co2 - pH_adm));
    alfa_IN = (pow(10, pK_a_IN - pH_adm))/(1.0 + pow(10, pK_a_IN - pH_adm));
    alfa_NH = 1.0/14000.0;  /* convert mgN/l into kmoleN/m3 */
    alfa_alk = -0.001;      /* convert moleHCO3/m3 into kmoleHCO3/m3 */
    alfa_NO = -1.0/14000.0; /* convert mgN/l into kmoleN/m3 */

	for (i = 0; i < 23; i++) {
     	utemp[i] = u[i];
        utemp2[i] = u[i];
    }
	
	for (i = 0; i < 32; i++)
		y[i] = 0.0; 

    /*================================================================================================*/
    /* Let CODdemand be the COD demand of available electron 
    * acceptors prior to the anaerobic digester, i.e. oxygen and nitrate */
   /* u[8]=u[8]+0.6*u[16]+u[17]*0.4+u[18]*0.2;/*MODIFIED*/
    y[134] = 0.086*(u[5]+u[20]+u[21]+u[4])+0.06*u[6]+u[9]+u[10]+u[11]+u[8]+u[16]+u[17]+u[18]+u[19];
    /*u[19]=u[19]+0.4*u[16]+u[17]*0.6+u[18]*0.8;/*ADDED*/
    u[5]=utemp[5]+utemp[20]+utemp[21];/*MODIFIED*/ /*AOB + NOB + AnAOB*/
    /*u[8]=u[8]+0.6*u[16]+0.4*u[18]+0.2*u[18];/*MODIFIED*/
    u[19]=u[19]+u[8]+u[16]+u[17]+u[18];
    CODdemand = u[7] + CODequiv_NO3*u[8]+CODequiv_NO2*u[16]+CODequiv_NO*u[17]+CODequiv_N2O*u[18];
    O2 = u[7];
    Y_H=mxGetPr(PAR)[29];
    n_Y=mxGetPr(PAR)[31];
    f_P=mxGetPr(PAR)[28];
    iN_XP = mxGetPr(PAR)[30];
    NO3 = u[8];
    NO = u[17];
    NO2 = u[16];
    N2O = u[18];
    utemp[11] = u[11] + u[5]*(iN_B-f_P*iN_XP);
    utemp[3] = u[3] + (1-f_P)*u[5];
    utemp[6] = u[6] + f_P*u[5];
    utemp[5] = 0;
    utemp[20] = 0 ;
    utemp[21] = 0;
    
        
   
     
    
    if ((u[8] = 0) && (u[7] = 0) && (u[16] = 0) && (u[17]=0) && (u[18]=0) && (u[19]=0)) {
         
         utemp[11] = utemp[11] + utemp[4]*(iN_B-f_P*iN_XP);
         utemp[3] =  utemp[3] + (1-f_P)*utemp[4];
         utemp[6] =  utemp[6] + f_P*utemp[4];
         utemp[4] = 0;
     
    }
    else {
        /*SS+ O2*/
    if ((u[1] < O2 /(1-Y_H)) && (u[9] < iN_B*Y_H/(1-Y_H)*O2)  ) {
       
        O2_consumedA = min( u[9]/(iN_B*Y_H/(1-Y_H)) , u[1]/(1/(1-Y_H)) );
          y[37] = 0.0;
    }
    
    else if((u[9] < iN_B*Y_H/(1-Y_H)*O2) && (u[1] >= O2 /(1-Y_H))  ) {
    
        O2_consumedA = u[9]/(iN_B*Y_H/(1-Y_H)); 
               y[37] = 1.0;  
    }
    
    else if ((u[1] < O2 /(1-Y_H)) && (u[9] >= iN_B*Y_H/(1-Y_H)*O2) ) {
    
        O2_consumedA = u[1]/(1/(1-Y_H));
          y[37] = 2.0;
    }
        
    else {
        
        O2_consumedA = O2;
            y[37] = 3.0;
    }
    
     O2remainingA = O2 - O2_consumedA ;
     utemp[1] = u[1] - O2_consumedA /(1-Y_H);
     NHremainingA = u[9] - iN_B*Y_H/(1-Y_H)*O2_consumedA;
     XBHremainingA = u[4] + Y_H/(1-Y_H)*O2_consumedA;
     
     y[33] = O2remainingA;
     y[34] = utemp[1];
     y[35] = NHremainingA;
     y[36] = XBHremainingA;
     y[38] = O2 /(1-Y_H);
     y[39] = iN_B*Y_H/(1-Y_H)*O2;
     y[40] = u[9] ;
     y[41] = u[1] ;
     y[42] = u[7] ;
  
     /*XS + O2*/
   if ((u[3] < O2remainingA /(1-Y_H)) && (u[9] < iN_B*Y_H/(1-Y_H)*O2remainingA) ) {
       
        O2_consumedB = min( u[9]/(iN_B*Y_H/(1-Y_H)) , u[3]/(1/(1-Y_H)) );
        y[53] = 4.0;
    }
    
    else if((u[9] < iN_B*Y_H/(1-Y_H)*O2remainingA) && (u[3] >= O2remainingA /(1-Y_H))  ) {
    
        O2_consumedB = u[9]/(iN_B*Y_H/(1-Y_H)); 
        y[53] = 5.0;  
    }
    
    else if ((u[3] < O2remainingA /(1-Y_H)) && (u[9] >= iN_B*Y_H/(1-Y_H)*O2remainingA)  ) {
    
        O2_consumedB = u[1]/(1/(1-Y_H));
        y[53] = 6.0;
    }
        
    else {
        
        O2_consumedB = O2remainingA;
        y[53] = 7.0;
    } 
     
        O2remainingB = O2remainingA - O2_consumedB ;
        utemp[3] = u[3] - O2_consumedB /(1-Y_H);
        NHremainingB = NHremainingA - iN_B*Y_H/(1-Y_H)*O2_consumedB;
        XBHremainingB = XBHremainingA + Y_H/(1-Y_H)*O2_consumedB;

     y[43] = O2remainingB;
     y[44] = utemp[3];
     y[45] = NHremainingB;
     y[46] = XBHremainingB;
     y[48] = O2remainingA /(1-Y_H);
     y[49] = iN_B*Y_H/(1-Y_H)*O2remainingA;
     y[50] = u[9] ;
     y[51] = u[3] ;
     y[52] = u[7] ;
     
          
     utemp[9]= NHremainingB;
     utemp[4] = XBHremainingB;
     
     
     /** PARALLEL NITROGEN OXIDE REDUCTION **/ 
     
     /*NO3 vs SS*/
     
     
     
     
     if ((utemp[1] < NO3*CODequiv_NO3 /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3)  ) {
       
        NO3_consumedA = min( utemp[9]/(iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[1]/(CODequiv_NO3/(1-n_Y*Y_H)) );
       y[58] = 8;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3) && (utemp[1] >= NO3*CODequiv_NO3 /(1-Y_H*n_Y)) ) {
    
        NO3_consumedA = utemp[9]/(iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[58] = 9;   
    }
    
     else if ((utemp[1] < NO3*CODequiv_NO3 /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3)  ) {
    
        NO3_consumedA = utemp[1]/(CODequiv_NO3/(1-n_Y*Y_H));
        y[58] = 10;
    }
        
    else {
        
        NO3_consumedA= NO3;
        y[58] = 11;
    }
      
        NO3remainingA = NO3 - NO3_consumedA ;
        utemp[1] = utemp[1] - NO3_consumedA *CODequiv_NO3/(1-n_Y*Y_H);
        NHremainingC = utemp[9] - iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3_consumedA;
        XBHremainingC = utemp[4] + CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3_consumedA;

     y[54] = NO3remainingA;
     y[55] = utemp[1];
     y[56] = NHremainingC;
     y[57] = XBHremainingC;
     y[59] = NO3remainingA /(1-n_Y*Y_H);
     y[60] = iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3remainingA;
     y[61] = u[9] ;
     y[62] = u[3] ;
     y[63] = NO3 ;
     
     utemp[9] = NHremainingC;
     utemp[4] = XBHremainingC;
          
      /*NO3 vs XS*/
     
      if ((utemp[3] < NO3remainingA*CODequiv_NO3 /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3remainingA)  ) {
       
        NO3_consumedB = min( utemp[9]/(iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[3]/(CODequiv_NO3/(1-n_Y*Y_H)) );
        y[68] = 12;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3remainingA) && (utemp[3] >= NO3remainingA*CODequiv_NO3 /(1-Y_H*n_Y)) ) {
    
        NO3_consumedB = utemp[9]/(iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[68] = 13;
     }
    
     else if ((utemp[3] < NO3remainingA*CODequiv_NO3 /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO3*Y_H*n_Y/(1-Y_H*n_Y)*NO3remainingA)  ) {
    
        NO3_consumedB = utemp[3]/(CODequiv_NO3/(1-n_Y*Y_H));
        y[68] = 14;
    }
        
    else {
        
        NO3_consumedB= NO3remainingA;
        y[68] = 15;
   }
      
     
     NO3remainingB = NO3remainingA - NO3_consumedB ;
     utemp[3] = utemp[3] - NO3_consumedB *CODequiv_NO3/(1-n_Y*Y_H);
     NHremainingD = utemp[9] - iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3_consumedB;
     XBHremainingD = utemp[4] + CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3_consumedB;
     
     y[64] = NO3remainingB;
     y[65] = utemp[3];
     y[66] = NHremainingC;
     y[67] = XBHremainingC;
     y[69] = NO3remainingB /(1-n_Y*Y_H);
     y[70] = iN_B*CODequiv_NO3*n_Y*Y_H/(1-n_Y*Y_H)*NO3remainingB;
     y[71] = u[9] ;
     y[72] = u[3] ;
     y[73] = NO3 ;
     
     utemp[9] = NHremainingD;
     utemp[4] = XBHremainingD;

     
     
     /*NO2 vs SS*/
     
     
     
     
     if ((utemp[1] < NO2*CODequiv_NO2 /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2)  ) {
       
        NO2_consumedA = min( utemp[9]/(iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[1]/(CODequiv_NO2/(1-n_Y*Y_H)) );
        y[78] = 16;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2) && (utemp[1] >= NO2*CODequiv_NO2 /(1-Y_H*n_Y))  ) {
    
        NO2_consumedA = utemp[9]/(iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[78] = 17;   
    }
    
     else if ((utemp[1] < NO2*CODequiv_NO2 /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2)  ) {
    
        NO2_consumedA = utemp[1]/(CODequiv_NO2/(1-n_Y*Y_H));
        y[78] = 18;
    }
        
    else {
        
        NO2_consumedA= NO2;
         y[78] = 19;
    }
      
        NO2remainingA = NO2 - NO2_consumedA ;
        utemp[1] = utemp[1] - NO2_consumedA *CODequiv_NO2/(1-n_Y*Y_H);
        NHremainingE = utemp[9] - iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2_consumedA;
        XBHremainingE = utemp[4] + CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2_consumedA;
     y[74] = NO2remainingA;
     y[75] = utemp[1];
     y[76] = NHremainingE;
     y[77] = XBHremainingE;
     y[79] = NO2remainingA /(1-n_Y*Y_H);
     y[80] = iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2remainingA;
     y[81] = u[9] ;
     y[82] = u[3] ;
     y[83] = u[16] ;
     
     utemp[9] = NHremainingE;
     utemp[4] = XBHremainingE;
          
      /*NO2 vs XS*/
     
      if ((utemp[3] < NO2remainingA*CODequiv_NO2 /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2remainingA)  ) {
       
        NO2_consumedB = min( utemp[9]/(iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[3]/(CODequiv_NO2/(1-n_Y*Y_H)) );
        y[88] = 20;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2remainingA) && (utemp[3] >= NO2remainingA*CODequiv_NO2 /(1-Y_H*n_Y))  ) {
    
        NO2_consumedB = utemp[9]/(iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[88] = 21;
     }
    
     else if ((utemp[3] < NO2remainingA*CODequiv_NO2 /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO2*Y_H*n_Y/(1-Y_H*n_Y)*NO2remainingA)  ) {
    
        NO2_consumedB = utemp[3]/(CODequiv_NO2/(1-n_Y*Y_H));
        y[88] =22;
    }
        
    else {
        
        NO2_consumedB= NO2remainingA;
        y[88] = 23;
   }
        
        NO2remainingB = NO2remainingA - NO2_consumedB ;
        utemp[3] = utemp[3] - NO2_consumedB *CODequiv_NO2/(1-n_Y*Y_H);
        NHremainingF = utemp[9] - iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2_consumedB;
        XBHremainingF = utemp[4] + CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2_consumedB;
 
     y[84] = NO2remainingB;
     y[85] = utemp[3];
     y[86] = NHremainingF;
     y[87] = XBHremainingF;
     y[89] = NO2remainingB /(1-n_Y*Y_H);
     y[90] = iN_B*CODequiv_NO2*n_Y*Y_H/(1-n_Y*Y_H)*NO2remainingB;
     y[91] = u[9] ;
     y[92] = u[3] ;
     y[93] = u[16] ;
     
     utemp[9] = NHremainingF;
     utemp[4] = XBHremainingF;
     
     
        /*NO vs SS*/
     
     
    
     
     if ((utemp[1] < NO*CODequiv_NO /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NO)  ) {
       
        NO_consumedA = min( utemp[9]/(iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[1]/(CODequiv_NO/(1-n_Y*Y_H)) );
        y[97] = 24;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NO) && (utemp[1] >= NO*CODequiv_NO /(1-Y_H*n_Y))  ) {
    
        NO_consumedA = utemp[9]/(iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[97] = 25;   
    }
    
     else if ((utemp[1] < NO*CODequiv_NO /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NO)  ) {
    
        NO_consumedA = utemp[1]/(CODequiv_NO/(1-n_Y*Y_H));
        y[97] = 26;
    }
        
    else {
        
        NO_consumedA= NO;
         y[97] = 27;
    }
      
        NOremainingA = NO - NO_consumedA ;
        utemp[1] = utemp[1] - NO_consumedA *CODequiv_NO/(1-n_Y*Y_H);
        NHremainingG = utemp[9] - iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NO_consumedA;
        XBHremainingG = utemp[4] + CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NO_consumedA;
        
        
     y[94] = NOremainingA;
     y[95] = NHremainingG;
     y[96] = XBHremainingG;
     y[98] = NOremainingA /(1-n_Y*Y_H);
     y[99] = iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NOremainingA;
     y[100] = u[9] ;
     y[101] = utemp[1] ;
     y[102] = u[17] ;
     
     utemp[9] = NHremainingG;
     utemp[4] = XBHremainingG;

     
         /*NO vs XS*/
     
     
     if ((utemp[3] < NOremainingA*CODequiv_NO /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NOremainingA)  ) {
       NO_consumedB = min( utemp[9]/(iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[3]/(CODequiv_NO/(1-n_Y*Y_H)) );
       
        y[107] = 28;
    }
    
     else if((utemp[9] < iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NOremainingA) && (utemp[3] >= NOremainingA*CODequiv_NO /(1-Y_H*n_Y))  ) {
    
        NO_consumedB = utemp[9]/(iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[107] = 29;   
    }
    
     else if ((utemp[3] < NOremainingA*CODequiv_NO /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_NO*Y_H*n_Y/(1-Y_H*n_Y)*NOremainingA)  ) {
    
        NO_consumedB = utemp[3]/(CODequiv_NO/(1-n_Y*Y_H));
        y[107] = 30;
    }
        
    else {
        
        NO_consumedB= NOremainingA;
         y[107] =31;
    }
      
        NOremainingB = NOremainingA - NO_consumedB ;
        utemp[3] = utemp[3] - NO_consumedB *CODequiv_NO/(1-n_Y*Y_H);
        NHremainingH = utemp[9] - iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NO_consumedB;
        XBHremainingH = utemp[4] + CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NO_consumedB;
        
        
     y[104] = NOremainingB;
     y[105] = NHremainingH;
     y[106] = XBHremainingH;
     y[108] = NOremainingB /(1-n_Y*Y_H);
     y[109] = iN_B*CODequiv_NO*n_Y*Y_H/(1-n_Y*Y_H)*NOremainingB;
     y[110] = u[9] ;
     y[111] = utemp[3] ;
     y[112] = u[17] ;
     
     utemp[9] = NHremainingH;
     utemp[4] = XBHremainingH;
    
     
     /*N2O + SS*/
     
      
     if ((utemp[1] < N2O*CODequiv_N2O /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2O)  ) {
       
        N2O_consumedA = min( utemp[9]/(iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[1]/(CODequiv_N2O/(1-n_Y*Y_H)) );
        y[116] = 32;
    }
    
     else if((utemp[9] < iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2O) && (utemp[1] >= N2O*CODequiv_N2O /(1-Y_H*n_Y))  ) {
    
        N2O_consumedA = utemp[9]/(iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[116] = 33;   
    }
    
     else if ((utemp[1] < N2O*CODequiv_N2O /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2O)  ) {
    
        N2O_consumedA = utemp[1]/(CODequiv_N2O/(1-n_Y*Y_H));
        y[116] = 34;
    }
        
    else {
        
        N2O_consumedA= N2O;
         y[116] = 35;
    }
      
        N2OremainingA = N2O - N2O_consumedA ;
        utemp[1] = utemp[1] - N2O_consumedA *CODequiv_N2O/(1-n_Y*Y_H);
        NHremainingI = utemp[9] - iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2O_consumedA;
        XBHremainingI = utemp[4] + CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2O_consumedA;
        
        
     y[113] = N2OremainingA;
     y[114] = NHremainingI;
     y[115] = XBHremainingI;
     y[117] = N2OremainingA /(1-n_Y*Y_H);
     y[118] = iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2OremainingA;
     y[119] = u[9] ;
     y[120] = utemp[1] ;
     y[121] = u[18] ;
     
     utemp[9] = NHremainingI;
     utemp[4] = XBHremainingI;
 
     
     
         /*N2O vs XS*/
     
     
     if ((utemp[3] < N2OremainingA*CODequiv_N2O /(1-Y_H*n_Y)) && (utemp[9] < iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2OremainingA)  ) {
      
        N2O_consumedB = min( utemp[9]/(iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)) , utemp[3]/(CODequiv_N2O/(1-n_Y*Y_H)) );
        y[125] = 36;
    }
    
     else if((utemp[9] < iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2OremainingA) && (utemp[3] >= N2OremainingA*CODequiv_N2O /(1-Y_H*n_Y))  ) {
    
        N2O_consumedB = utemp[9]/(iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)); 
        y[125] = 37;   
    }
    
     else if ((utemp[3] < N2OremainingA*CODequiv_N2O /(1-Y_H*n_Y)) && (utemp[9] >= iN_B*CODequiv_N2O*Y_H*n_Y/(1-Y_H*n_Y)*N2OremainingA)  ) {
    
        N2O_consumedB = utemp[3]/(CODequiv_N2O/(1-n_Y*Y_H));
        y[125] = 38;
    }
        
    else {
        
        N2O_consumedB= N2OremainingA;
         y[125] =39;
    }
      
        N2OremainingB = N2OremainingA - N2O_consumedB ;
        utemp[3] = utemp[3] - N2O_consumedB *CODequiv_N2O/(1-n_Y*Y_H);
        NHremainingL = utemp[9] - iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2O_consumedB;
        XBHremainingL = utemp[4] + CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2O_consumedB;
        
        
     y[122] = N2OremainingB;
     y[123] = NHremainingL;
     y[124] = XBHremainingL;
     y[126] = N2OremainingB /(1-n_Y*Y_H);
     y[127] = iN_B*CODequiv_N2O*n_Y*Y_H/(1-n_Y*Y_H)*N2OremainingB;
     y[128] = u[9] ;
     y[129] = utemp[3] ;
     y[130] = u[18] ;
     
     utemp[9] = NHremainingL;
     utemp[4] = XBHremainingL;
 
    }
     
    
    
         /*XBH DECAY*/
         utemp[11] = utemp[11] + utemp[4]*(iN_B-f_P*iN_XP);
         utemp[3] =  utemp[3] + (1-f_P)*utemp[4];
         utemp[6] =  utemp[6] + f_P*utemp[4];
         utemp[4] = 0;
     
         y[131]= utemp[11];
         y[132]= utemp[3];
         y[133] = utemp[6];
     
         u[11] = utemp[11];
         u[10] = utemp[10];
         
         y[135] =  0.086*(utemp[5]+utemp[20]+utemp[21]+utemp[4])+0.06*utemp[6]+utemp[9]+utemp[10]+utemp[11]+utemp[8]+utemp[16]+utemp[17]+utemp[18]+u[19];
    
        /*if ((u[3] <= O2remainingA /(1-Y_H)) || (NHremainingA <= iN_B*Y_H/(1-Y_H)*O2remainingA)) {
            |
        O2_consumedB = min( NHremainingA/(iN_B*Y_H/(1-Y_H)) , u[3]/(1/(1-Y_H)) );
        O2remainingB = O2remainingA - O2_consumedB ;
        utemp[3] = u[3] - O2_consumedB /(1-Y_H);
        NHremainingB = NHremainingA - iN_B*Y_H/(1-Y_H)*O2_consumedB;
        XBHremainingB = XBHremainingA + Y_H/(1-Y_H)*O2_consumedB;
           
        }
    }
        
    /*CODdemand_remaining1 = CODdemand - (1 - Y_H)*u[1] ; 
    utemp[1] = 0;
    deltaBIOMASS1 = Y_H * u[1];
    deltaNH1 = - deltaBIOMASS1 * iN_B;
    
    CODdemand_remaining2 = CODdemand_remaining1 - (1 - Y_H)*u[3];
    utemp[3] = 0;
    deltaBIOMASS2 = Y_H*u[3]; 
    deltaNH2 =  - deltaBIOMASS2 * iN_B;
    
    NtotIN = u[9] + u[10] + u[11];
    Nremaining = NtotIN +(deltaNH1 + deltaNH2);
        
    y[33] = Nremaining;
    y[34] = deltaNH1;
    y[35] = deltaNH2;
    
    CODdemand_remaining3 = CODdemand_remaining2 - (u[4] + deltaBIOMASS1 + deltaBIOMASS2);
    u[4] = 0;
    deltaBIOMASS3 = -(u[4] + deltaBIOMASS1 + deltaBIOMASS2);
    deltaNH3 = -deltaBIOMASS3*iN_B;
    
    CODdemand_remaining4 = CODdemand_remaining3 - u[5];
    u[5] = 0 ;
    deltaBIOMASS3bis = -u[5];
    deltaNH4 = -deltaBIOMASS3bis * iN_B;
    
    deltaNHtot = deltaNH1 + deltaNH2 + deltaNH3 + deltaNH4;
    
    y[36] = deltaNHtot;
    y[37] = u[9] + u[10] + u[11] - y[36];
    
    utemp[7] = 0;
    utemp[8] = 0;
    utemp[16]= 0;
    utemp[17]=0;
    utemp[18] =0;

    /* if extreme detail was used then some extra NH4 would be transformed
    * into N bound in biomass and some biomass would be formed when
    * removing the CODdemand (based on the yield). But on a total COD balance 
	* approach the below is correct (neglecting the N need for biomass growth)
    * The COD is reduced in a hierarchical approach in the order: 
    * 1) SS; 2) XS; 3) XBH; 4) XBA. It is no real improvement to remove SS and add
    * biomass. The net result is the same.*/
	
	/*if (CODdemand > u[1]) {	/* check if COD demand can be fulfilled by SS*/
  	/*	remaina = CODdemand - u[1];
  		utemp[1] = 0.0;
  		if (remaina > u[3]) {	/* check if COD demand can be fulfilled by XS*/
    	/*	remainb = remaina - u[3];
    		utemp[3] = 0.0;
    		if (remainb > u[4]) {	/* check if COD demand can be fulfilled by XBH */
      		/*	remainc = remainb - u[4];
                utemp[9] = utemp[9] + u[4]*iN_B;
      			utemp[4] = 0.0;
      			if (remainc > u[5]) {	/* check if COD demand can be fulfilled by XBA */
        	/*		remaind = remainc - u[5];
                    utemp[9] = utemp[9] + u[5]*iN_B;
        			utemp[5] = 0.0;
                    utemp[7] = remaind;
					/* if here we are in trouble, carbon shortage: an error printout should be given */
                    /* and execution stopped */
      		/*	}
        		else {		/* reduced all COD demand by use of SS, XS, XBH and XBA */
       /* 			utemp[5] = u[5] - remainc;
                    utemp[9] = utemp[9] + remainc*iN_B;
                }
            }
    		else {		/* reduced all COD demand by use of SS, XS and XBH */
      /*			utemp[4] = u[4] - remainb;
                utemp[9] = utemp[9] + remainb*iN_B;
            }
		}
  		else {		/* reduced all COD demand by use of SS and XS */
  /*  		utemp[3] = u[3] - remaina; 
		}
	}
	else {		/* reduced all COD demand by use of SS */
/*  		utemp[1] = u[1] - CODdemand;
	}

    /*================================================================================================*/
	/* SS becomes part of amino acids when transformed into ADM
    * and any remaining SS is mapped to monosacharides (no N contents)
    * Enough SND must be available for mapping to amino acids */
    
	sorgn = u[10]/fnaa;     /* Saa COD equivalent to SND */

	if (sorgn >= utemp[1]) {	/* not all SND-N in terms of COD fits into amino acids */
  		y[1] = utemp[1];        /* map all SS COD into Saa */
  		utemp[10] = utemp[10] - utemp[1]*fnaa;	/* excess SND */
  		utemp[1] = 0.0;			/* all SS used */
	}
	else {                      /* all SND-N fits into amino acids */
 	 	y[1] = sorgn;           /* map all SND related COD into Saa */
  		utemp[1] = utemp[1] - sorgn;		/* excess SS, which will become sugar in ADM1 i.e. no nitrogen association */
  		utemp[10] = 0.0;		/* all SND used */
	}

    /*================================================================================================*/
	/* XS becomes part of Xpr (proteins) when transformed into ADM
    * and any remaining XS is mapped to Xch and Xli (no N contents)
    * Enough XND must be available for mapping to Xpr */
    
	xorgn = u[11]/fnaa;     /* Xpr COD equivalent to XND */

	if (xorgn >= utemp[3]) {        /* not all XND-N in terms of COD fits into Xpr */
  		xprtemp = utemp[3];         /* map all XS COD into Spr */
  		utemp[11] = utemp[11] - utemp[3]*fnaa;	/* excess XND */
  		utemp[3] = 0.0;             /* all XS used */
        xlitemp = 0.0;
        xchtemp = 0.0;
	}
	else {                      /* all XND-N fits into Xpr */
 	 	xprtemp = xorgn;        /* map all XND related COD into Xpr */
        xlitemp = frlixs*(utemp[3] - xorgn);    /* part of XS COD not associated with N */
        xchtemp = (1.0 - frlixs)*(utemp[3] - xorgn);    /* part of XS COD not associated with N */
        utemp[3] = 0.0;           /* all XS used */
  		utemp[11] = 0.0;		/* all XND used */
	}

    /*================================================================================================*/
    /* Biomass becomes part of Xpr and XI when transformed into ADM
	* and any remaining XBH and XBA is mapped to Xch and Xli (no N contents)
	* Remaining XND-N can be used as nitrogen source to form Xpr */
    
    biomass = utemp[4] + utemp[5]+utemp[20]+utemp[21]; /*MODIFIED*/
    biomass_nobio = biomass*(1.0 - frxs_adm);   /* part which is mapped to XI */
    biomass_bioN = (biomass*iN_B - biomass_nobio*fxni);
    if (biomass_bioN < 0.0) {
        /* Problems: if here we should print 'ERROR: not enough biomass N to map the requested inert part' */
        /* and execution stopped */
    }
    if ((biomass_bioN/fnaa) <= (biomass - biomass_nobio)) {
        xprtemp2 = biomass_bioN/fnaa;   /* all biomass N used */
        remainCOD = biomass - biomass_nobio - xprtemp2;
        if ((utemp[11]/fnaa) > remainCOD) {  /* use part of remaining XND-N to form proteins */
            xprtemp2 = xprtemp2 + remainCOD;
            utemp[11] = utemp[11] - remainCOD*fnaa;
            remainCOD = 0.0;
            utemp[4] = 0.0;
            utemp[5] = 0.0;
            utemp[20]=0.0;
            utemp[21]=0.0;
        }
        else {       /* use all remaining XND-N to form proteins */
            xprtemp2 = xprtemp2 + utemp[11]/fnaa;
            remainCOD = remainCOD - utemp[11]/fnaa;
            utemp[11] = 0.0;
        }
        xlitemp2 = frlibac*remainCOD;        /* part of the COD not associated with N */
        xchtemp2 = (1.0 - frlibac)*remainCOD;  /* part of the COD not associated with N */
    }
    else {
        xprtemp2 = biomass - biomass_nobio; /* all biomass COD used */
        utemp[11] = utemp[11] + biomass*iN_B - biomass_nobio*fxni - xprtemp2*fnaa; /* any remaining N in XND */
    }
    utemp[4] = 0.0;
    utemp[5] = 0.0;
    utemp[20] = 0.0; /*added*/
    utemp[21]=0.0;
    /*================================================================================================*/
    /* direct mapping of XI and XP to ADM1 XI (if fdegrade_ad = 0)
	* assumption: same N content in both ASM1 and ADM1 particulate inerts */
    
    inertX = (1-fdegrade_adm)*(utemp[2] + utemp[6]);

    /* special case: IF part of XI and XP in the ASM can be degraded in the AD
    * we have no knowledge about the contents so we put it in as composits (Xc)
	* we need to keep track of the associated nitrogen
	* N content which may be different, take first from XI&XP-N, then XND-N, then SND-N,
	* then SNH. A similar principle could be used for other states. */

    xc = 0.0;
    xlitemp3 = 0.0;
    xchtemp3 = 0.0;
    if (fdegrade_adm > 0)  {
        noninertX = fdegrade_adm*(utemp[2] + utemp[6]);
        if (fxni < fnxc) {   /* N in XI&XP(ASM) not enough */
            xc = noninertX*fxni/fnxc;
            noninertX = noninertX - noninertX*fxni/fnxc;
            if (utemp[11] < (noninertX*fnxc)) {   /* N in XND not enough */
                xc = xc + utemp[11]/fnxc;
                noninertX = noninertX - utemp[11]/fnxc;
                utemp[11] = 0.0;
                if (utemp[10] < (noninertX*fnxc)) {   /* N in SND not enough */
                    xc = xc + utemp[10]/fnxc;
                    noninertX = noninertX - utemp[10]/fnxc;
                    utemp[10] = 0.0; 
                    if (utemp[9] < (noninertX*fnxc)) {   /* N in SNH not enough */
                        xc = xc + utemp[9]/fnxc;
                        noninertX = noninertX - utemp[9]/fnxc;
                        utemp[9] = 0.0;
                        /* Should be a WARNING printout: Nitrogen shortage when converting biodegradable XI&XP 
                        * Putting remaining XI&XP as lipids (50%) and carbohydrates (50%) */ 
                        xlitemp3 = 0.5*noninertX;
                        xchtemp3 = 0.5*noninertX;
                        noninertX = 0.0;
                        }
                    else {   /* N in SNH enough for mapping */
                        xc = xc + noninertX;
                        utemp[9] = utemp[9] - noninertX*fnxc;
                        noninertX = 0.0;
                        }
                    }
                else  {   /* N in SND enough for mapping */
                    xc = xc + noninertX;
                    utemp[10] = utemp[10] - noninertX*fnxc;
                    noninertX = 0.0;
                    }
                }
            else  {   /* N in XND enough for mapping */
                xc = xc + noninertX;
                utemp[11] = utemp[11] - noninertX*fnxc;
                noninertX = 0.0;
                }
            }
        else  {   /* N in XI&XP(ASM) enough for mapping */
            xc = xc + noninertX;
            utemp[11] = utemp[11] + noninertX*(fxni-fnxc);   /* put remaining N as XND */
            noninertX = 0;
            }
    }

    /*================================================================================================*/
    /* Mapping of ASM SI to ADM1 SI
	* N content may be different, take first from SI-N, then SND-N, then XND-N,
	* then SNH. Similar principle could be used for other states. */
    
    inertS = 0.0;
    if (fsni < fsni_adm) {   /* N in SI(ASM) not enough */
        inertS = utemp[0]*fsni/fsni_adm;
        utemp[0] = utemp[0] - utemp[0]*fsni/fsni_adm;
        if (utemp[10] < (utemp[0]*fsni_adm)) {    /* N in SND not enough */
            inertS = inertS + utemp[10]/fsni_adm;
            utemp[0] = utemp[0] - utemp[10]/fsni_adm;
            utemp[10] = 0.0;
            if (utemp[11] < (utemp[0]*fsni_adm)) {   /* N in XND not enough */
                inertS = inertS + utemp[11]/fsni_adm;
                utemp[0] = utemp[0] - utemp[11]/fsni_adm;
                utemp[11] = 0.0; 
                if (utemp[9] < (utemp[0]*fsni_adm)) {   /* N in SNH not enough */
                    inertS = inertS + utemp[9]/fsni_adm;
                    utemp[0] = utemp[0] - utemp[9]/fsni_adm;
                    utemp[9] = 0.0;
                    /* Here there shpuld be a warning printout: Nitrogen shortage when converting SI 
                    * Putting remaining SI as monosacharides */
                    utemp[1] = utemp[1] + utemp[0];
                    utemp[0] = 0.0;
                    }
                else  {   /* N in SNH enough for mapping */
                    inertS = inertS + utemp[0];
                    utemp[9] = utemp[9] - utemp[0]*fsni_adm;
                    utemp[0] = 0.0; 
                    }
                }
            else  {   /* N in XND enough for mapping */
                inertS = inertS + utemp[0];
                utemp[11] = utemp[11] - utemp[0]*fsni_adm;
                utemp[0] = 0.0;
                }
            }
        else  {  /* N in SND enough for mapping */
            inertS = inertS + utemp[0];
            utemp[10] = utemp[10] - utemp[0]*fsni_adm;
            utemp[0] = 0.0;
            }
        }
    else  {   /* N in SI(ASM) enough for mapping */
        inertS = inertS + utemp[0];
        utemp[10] = utemp[10] + utemp[0]*(fsni-fsni_adm);   /* put remaining N as SND */
        utemp[0] = 0.0;
        }

    /*================================================================================================*/
    /* Define the outputs including charge balance */
    
    y[0] = utemp[1]/1000.0;/*Monosaccarides*/
    y[1] = y[1]/1000.0;
    y[10] = (utemp[9] + utemp[10] + utemp[11])/14000.0;
    y[11] = inertS/1000.0;
    y[12] = xc/1000.0;
    y[13] = (xchtemp + xchtemp2 + xchtemp3)/1000.0;
    y[14] = (xprtemp + xprtemp2)/1000.0;
    y[15] = (xlitemp + xlitemp2 + xlitemp3)/1000.0;
    y[23] = (biomass_nobio + inertX)/1000.0;
    y[26] = u[14];  /* flow rate */
    y[27] = T_op - 273.15;   /* temperature, degC */
    /*y[28] = u[16];  /* dummy state */
    y[28] = 0;
    /*y[29] = u[17];  /* dummy state */
    y[29] = 0;
    /*y[30] = u[18];  /* dummy state */
    y[30] = 0;
    y[31] = u[19];  /* dummy state */
    y[32] = 0;  /* dummy state */

    /* charge balance, output S_IC */
    y[9] = (((utemp2[8]+utemp2[16])*alfa_NO + utemp2[9]*alfa_NH + utemp2[12]*alfa_alk) - (y[3]*alfa_va + y[4]*alfa_bu + y[5]*alfa_pro + y[6]*alfa_ac + y[10]*alfa_IN))/alfa_co2;/*MODIFIED*/

    /* calculate anions and cations based on full charge balance including H+ and OH- */
    ScatminusSan = y[3]*alfa_va + y[4]*alfa_bu + y[5]*alfa_pro + y[6]*alfa_ac + y[10]*alfa_IN + y[9]*alfa_co2 + pow(10, (-pK_w + pH_adm)) - pow(10, -pH_adm);
    
    if (ScatminusSan > 0)  {
        y[24] = ScatminusSan;
        y[25] = 0.0;
        }
    else  {
        y[24] = 0.0;
        y[25] = -1.0*ScatminusSan;
        }

    /* Finally there should be a input-output mass balance check here of COD and N */
    
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

