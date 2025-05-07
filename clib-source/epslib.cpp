#include <iostream> 
#include <cstdio>
#include <complex.h>
#include <math.h>
#include <quadmath.h>
#include <string>
#include <cstdlib>
#include <float.h>
#include  <cmath>
#include <fstream>
#include <cerrno>
#include <cfenv>
// #pragma STDC FENV_ACCESS ON
#  include "Faddeeva.hh"

using namespace std;
typedef complex< double> dcomp;
//typedef complex<long  double> ldcomp;
#define MAXOSC  250// max number of oscillators, should match what is happening in calculate.py
#define MAXGOS  10 // Max number of core levels described by GOS, should match what is happening in calculate.py
#define MAXBELKACEM 5// maximum number of Belkacem oscillators, should match what is happening in calculate.py
#define MAXKANEKO 250 // Max number of core levels described by Kaneko, should match what is happening in calculate.py

const double pi = 3.141592653589793238463;
const double Hartree = 27.211;
const double BohrRadius  = 0.529177;
const double C = 137.036;
const double MassProton =1836.15;
double Projectile_Mass;
bool proton;  // true = proton else electron
bool MottCorrection;
bool debugmessageprinted=false;
bool DebugMode=false;

double Ai[MAXOSC], gammai[MAXOSC], wi[MAXOSC], alphai[MAXOSC],gapi[MAXOSC];
double AiGOS[MAXGOS ],Edgei_GOS[MAXGOS],Zi_GOS[MAXGOS],GOS_ScalingFactor[MAXGOS][80];
double TaucNormalisationq0[MAXOSC];
double maxEnergyDensityEffect;
const double GOSScalingStepSize=0.2; // used in GOS_Scaling_init
double Ai_Belkacem[MAXBELKACEM],wi_Belkacem[MAXBELKACEM],gammai_Belkacem[MAXBELKACEM];
bool  Full_dispersion,Local_Desep, Dispersion_relativistic, kk_variable_stepsize;
int n_i_GOS[MAXGOS ], l_i_GOS[MAXGOS ];

double N_Kaneko[MAXKANEKO ],Edge_ArchubiKaneko[MAXKANEKO],Q_Kaneko[MAXKANEKO],Width_Kaneko[MAXKANEKO];
double w_pl_l[MAXKANEKO],w_pl_0[MAXKANEKO],A_l_Kaneko[MAXKANEKO],A_0_Kaneko[MAXKANEKO],gamma_Kaneko[MAXKANEKO];
int l_Kaneko[MAXKANEKO];

double UnitCellDensity, PL_GOS;
double epsbkg;
double w_global, theta_global; 
int I_Os_global;
double quanc8result; double quanc8errest; int quanc8nofun; double quanc8flag;  // not sure why I can not pass addresses here and have ot make the output a global variable
double E_0,p_0,v_0, rel_cor_factor,Egerton_rel_cor_factor, gamma_rel, beta_r, b_zero;

double rest_mass_energy;

bool modelMerminLL,Apply_Mermin_Correction, modelDrude, modelDL,modelTaucLorentz, DirectMethod,  ApplySumRuleToGOS, OriginalKaneko, AddELF;
int ExchangeCorrection;
double abserr, relerr;
double FirstEnergy ,StepSize_eV,StepSize,Stepsize_qplot ,LastMomentum; 
double q_lower,q_upper;
int NStep;
double lin_cont_deltaE; // for the integration of the DIIMFP for stopping, also Tauc normalisation
double q2maxFactor, q2surfmaxFactor;
double q_transition, BE_for_exchange;
double theta_max;

double mom_limit_lower[30], mom_limit_upper[30]; 
double FSumRule, BetheSumRule,GOSBetheSumRule,KKSumRule,A_ScalingFactor,w_ScalingFactor;
double theta0,theta1,sigma,PIcoef1,PIcoef2,PIcoef3,surf_ex_factor, fraction_DIIMFP;  // for reels spectrum
dcomp  oneoverDL(double q, double current_w, double  gamma, double w, double alpha, double U);  // function declaration we generally use function after definition.  this is the exception, hence needs a declaration
dcomp  oneoverTL(double q, double current_w, double  C, double E0, double alpha, double Egap);


int check_for_exceptions()
 {
    int fe;
    /* testing multiple exceptions: */
    fe = fetestexcept (FE_ALL_EXCEPT);
    if ((fe !=0)&&(fe != 32) )
    {
          if (fe & FE_DIVBYZERO) printf("FE_DIVBYZERO\n");
          if (fe & FE_INVALID)   printf("FE_INVALID\n");
          if (fe & FE_OVERFLOW)  printf("FE_OVERFLOW\n");
          if (fe & FE_UNDERFLOW) printf("FE_UNDERFLOW\n");
          return fe;
  }
 return 0;
}
void my_perror(const char* str)
{
     if (DebugMode)
     {
         if(errno !=0) perror(str);
         check_for_exceptions();
         errno=0;
    }
}

double Poisson(int  n_occurence, double Plambda)
{  double result;
    if (DebugMode) my_perror("before lgamma\n");
    result= exp(n_occurence*log(Plambda)-Plambda-lgamma(n_occurence+1.0));
     if (DebugMode) my_perror("after lgamma\n");
    return result;

}

inline double recoil_energy_free_electron( double q) // calculates the energy a free electron gets after it absorbes momentum q
{
    double Q_recoil;
     if(q < 3.0)  // nonrelativistic is good enough here, relativistic formula  not stable at small q
    {
        Q_recoil = q*q/2.0;
    }
    else
    {
        Q_recoil= sqrt(C * C * q * q + C * C * C * C) - C * C;
    }
    return Q_recoil;
}

inline dcomp oneoverchi_from_chi(dcomp chi)    //works both ways also gets chi from oneoverchi
{ 
    dcomp eps,oneovereps, oneoverchi;
    eps=chi+dcomp(1.0,0.0);
    oneovereps=dcomp(1.0,0)/eps;
    oneoverchi= oneovereps-dcomp(1.0,0.0);
    return oneoverchi;
}

double  GOSx(int n, int  l, const double Z, double  dE_min, double q, double  w) //results GOS PER ELECTRON
																		   // Z = 0 uses dE_min as binding energy
{           

	double Zs, ne, Wl, Ql, Ql2, Ql3, A, bl, dummy, kH, kH2, Qrecoil;
	double c[10]; 
	int j, jmax;
    if (w < dE_min) return 0.0;
	Zs = 0.0; 
	if ((n == 1)&&(Z==1))  Zs = Z;
    else if(n == 1) Zs = Z - 0.3;
	else if (n == 2)  Zs = Z - 4.15;
	else if (n == 3)
	{
		if ((l == 0) || (l == 1))  Zs = Z - 11.25;
		if (l == 2) Zs = Z - 21.15;
	}
	if (Z == 0.0)  Zs = n*sqrt(dE_min / 0.5);

	if (Zs <= 0.0) return 0.0;
	ne = 2.0 * (2.0 * l + 1.0);
 

	//Ql = (q / Zs)*(q / Zs);  original
    if(Dispersion_relativistic)
    {
        Qrecoil = recoil_energy_free_electron(q);
    }
    else
    {
        Qrecoil = q * q / 2.0;
    }
    Ql = 2.0 * Qrecoil / ( Zs * Zs );
	Wl = w / (0.5*Zs*Zs);
	kH2 = Wl - 1.0 / (n*n);

	if (kH2 > 0.0)
	{  
		kH = sqrt(kH2);
		bl = atan(2.0 * kH / n / (Ql - Wl + 2.0 / (n*n)));
       
		if (bl < 0)  bl = bl + pi;
        double tmp1=0.0;
        if(kH > 0.02) tmp1= exp(-2.0 * pi / kH);// in this way no underflow errors
		A = 16.0 * pow(2.0 / n, 3.0)*exp(-2.0 / kH*bl) / (1.0 - tmp1) / pow((Ql - Wl)*(Ql - Wl) + 4.0 / (n*n)*Ql, 2.0 * n + 1.0);
      
	}
	else
	{ 
		kH = sqrt(-kH2);
		double tmp = (Ql - Wl + 2.0 / (n*n) + 2.0 * kH / n) / (Ql - Wl + 2.0 / (n*n) - 2.0 * kH / n);
		bl = -1.0 / kH*log(tmp);
        
		A = 16.0 * pow(2.0 / n, 3.0)*exp(bl) / pow((Ql - Wl)*(Ql - Wl) + 4.0 / (n* n)*Ql, 2.0 * n + 1.0);
	}

	dummy = 0.0;
	jmax = -1;
  
	if (n == 1)
	{
		c[0] = (Ql + Wl / 3);
		jmax = 0;
	}
	if ((n == 2) && (l == 0))
	{
		Ql2 = Ql*Ql;
		c[0] = (19.0 / 60.0 + 8.0 / 15.0 * Ql)*Ql2;
		c[1] = 1.0 / 15.0 * (Ql + 1)*Ql;
		c[2] = 1.0 / 30.0 * (1.0 - 40.0 * Ql)*Ql;
		c[3] = 2.0 / 3.0 * Ql;
		c[4] = 1.0 / 4.0 + 4.0 / 3.0 * Ql;
		c[5] = 1.0 / 3.0;
		jmax = 5;
	}

	if ((n == 2) && (l == 1))
	{
		Ql2 = Ql*Ql;
		c[0] = (17.0 / 20.0 + 4.0 / 5.0 * Ql)*Ql2; // according to Spannish appendix correcting misprint Serra
												   //	c[0] = (17.0 / 12.0 + 4.0 / 5.0 * Ql)*Ql2; // according to Sera
		c[1] = (1.0 / 10.0 + 34.0 / 15.0 * Ql)*Ql;
		c[2] = 1.0 / 30.0 * (49.0 + 120.0 * Ql)*Ql;
		c[3] = 1.0 / 6.0 + 2.0 * Ql;
		c[4] = 1.0 / 4.0;
		jmax = 4;
            
	}

	if ((n == 3) && (l == 0))
	{
		Ql2 = Ql*Ql;
		Ql3 = Ql2*Ql;
		c[0] = (528384.0 / 502211745.0 + 561152.0 / 55801305.0 * Ql + 68608.0 / 6200145.0 * Ql2)*Ql3;
		c[1] = (32768.0 / 167403915.0 + 2048.0 / 413343.0 * Ql - 7424.0 / 6200145.0 * Ql2)*Ql2;
		c[2] = (8192.0 / 6200145.0 - 63488.0 / 1240029.0 * Ql - 17408.0 / 98415.0 * Ql2)*Ql2;
		c[3] = (256768.0 / 6200145.0 + 256.0 / 729.0 * Ql)*Ql2;
		c[4] = (11008.0 / 885735.0 + 30976.0 / 98415.0 * Ql + 15744.0 / 10935.0 * Ql2)*Ql;
		c[5] = (2816.0 / 32805.0 - 10912.0 / 10935.0 * Ql)*Ql;
		c[6] = (128.0 / 19683.0 + 1024.0 / 10935.0 * Ql - 64.0 / 27.0 * Ql2);
		c[7] = (208.0 / 2187.0 + 80.0 / 81.0 * Ql);
		c[8] = (32.0 / 81.0 + 4.0 / 3.0 * Ql);
		c[9] = 1.0 / 3.0;
		jmax = 9;
	}


	if ((n == 3) && (l == 1))
	{
		Ql2 = Ql*Ql;
		Ql3 = Ql2*Ql;
		c[0] = (495616.0 / 167403915.0 + 443392.0 / 18600435.0 * Ql + 8192.0 / 413343.0 * Ql2)*Ql3;
		c[1] = (65536.0 / 167403915.0 + 546304.0 / 18600435.0 * Ql + 38912.0 / 413343.0 * Ql2)*Ql2;
		c[2] = (274944.0 / 18600435.0 + 135424.0 / 2066715.0 * Ql + 2048.0 / 6561.0 * Ql2)*Ql2;
		c[3] = (4096.0 / 2657205.0 + 15872.0 / 137781.0 * Ql - 22528.0 / 32805.0 * Ql2)*Ql;
		c[4] = (512.0 / 32805.0 + 2368.0 / 10935.0 * Ql - 512.0 / 243.0 * Ql2)*Ql;
		c[5] = (8992.0 / 32805.0 + 1664.0 / 729.0 * Ql)*Ql;
		c[6] = (224.0 / 6561.0 + 560.0 / 243.0 * Ql + 128.0 / 27.0 * Ql2);
		c[7] = (208.0 / 729.0 + 64.0 / 27.0 * Ql);
		c[8] = 8.0 / 27.0;
		jmax = 8;
	}


	if ((n == 3) && (l == 2))
	{
		Ql2 = Ql*Ql;
		Ql3 = Ql2*Ql;
		c[0] = (253952.0 / 55801305.0 + 904192.0 / 55801305.0 * Ql + 131072.0 / 6200145.0 * Ql2)*Ql3;
		c[1] = (45056.0 / 167403915.0 + 1440256.0 / 18600435.0 * Ql + 149504.0 / 6200145.0 * Ql2)*Ql2;
		c[2] = (140800.0 / 3720087.0 + 4211968.0 / 6200145.0 * Ql + 32768.0 / 98415.0 * Ql2)*Ql2;
		c[3] = (2048.0 / 885735.0 + 2657792.0 / 6200145.0 * Ql + 48128.0 / 32805.0 * Ql2)*Ql;
		c[4] = (75008.0 / 885735.0 + 193472.0 / 98415.0 * Ql + 8192.0 / 3645.0 * Ql2)*Ql;
		c[5] = (256.0 / 59049.0 + 6304.0 / 10935.0 * Ql + 22912.0 / 10935.0 * Ql2);
		c[6] = (736.0 / 19683.0 + 6416.0 / 10935.0 * Ql);
		c[7] = 80.0 / 2187.0;
		jmax = 7;
	}

	for (j = 0; j < (jmax + 1); j++)
	{
		dummy = dummy + c[j] * pow(Wl - Ql, j);
	}

	double result = A*dummy * 2.0 * w / (0.5*(Zs*Zs)*(0.5*(Zs*Zs))) / ne;
	return result;
}


void quanc8(double(*fun)(double), double a, double b, double abserr, double relerr)
/*
estimate the integral of fun(x) from a to b to a user provided tolerance.
an automatic adaptive routine based on the 8-panel newton-cotes rule.

input:
fun     the name of the integrand function subprogram fun(x).
a       the lower limit of integration.
b       the upper limit of integration.(b may be less than a.)
relerr  a relative error tolerance. (should be non-negative)
abserr  an absolute error tolerance. (should be non-negative)

output:
result  an approximation to the integral hopefully satisfying the
least stringent of the two error tolerances.
quanc8errest  an estimate of the magnitude of the actual error.
quanc8nofun   the number of function values used in calculation of result.
quanc8flag    a reliability indicator.  if quanc8flag is zero, then result
probably satisfies the error tolerance.  if quanc8flag is
xxx.yyy , then  xxx = the number of intervals which have
not converged and  0.yyy = the fraction of the interval
left to do when the limit on  quanc8nofun  was approached.

comments:
Alex Godunov (February 2007)
the program is based on a fortran version of program quanc8.f
*/
{
	double w0, w1, w2, w3, w4, area, x0, f0, stone, cor11, temp;
	double qprev,  tolerr;
	double qright[32], f[17], x[17], fsave[9][31], xsave[9][31];
	// double dabs,dmax1;
	int    levmin, levmax, levout, nomax, nofin, lev, nim, i, j;
	int    key;

	//  ***   stage 1 ***   general initialization

	levmin = 1;
	levmax = 30; // was 30
	levout = 6;
    //nomax = 5000;
	nomax = 250000;  // was 5000
	nofin = nomax - 8 * (levmax - levout + 128);
	//  trouble when quanc8nofun reaches nofin

	w0 = 3956.0 / 14175.0;
	w1 = 23552.0 / 14175.0;
	w2 = -3712.0 / 14175.0;
	w3 = 41984.0 / 14175.0;
	w4 = -18160.0 / 14175.0;

	//  initialize running sums to zero.

	quanc8flag = 0.0;
	quanc8result = 0.0;
	cor11 = 0.0;
	quanc8errest = 0.0;
	area = 0.0;
	quanc8nofun = 0;
	if (a == b) return;

	//  ***   stage 2 ***   initialization for first interval

	lev = 0;
	nim = 1;
	x0 = a;
	x[16] = b;
	qprev = 0.0;
	f0 = fun(x0);
	stone = (b - a) / 16.0;
	x[8] = (x0 + x[16]) / 2.0;
	x[4] = (x0 + x[8]) / 2.0;
	x[12] = (x[8] + x[16]) / 2.0;
	x[2] = (x0 + x[4]) / 2.0;
	x[6] = (x[4] + x[8]) / 2.0;
	x[10] = (x[8] + x[12]) / 2.0;
	x[14] = (x[12] + x[16]) / 2.0;
	for (j = 2; j <= 16; j = j + 2)
	{
		f[j] = fun(x[j]);
	}
	quanc8nofun = 9;

	//  ***   stage 3 ***   central calculation

	while (quanc8nofun <= nomax)
	{
		x[1] = (x0 + x[2]) / 2.0;
		f[1] = fun(x[1]);
		for (j = 3; j <= 15; j = j + 2)
		{
			x[j] = (x[j - 1] + x[j + 1]) / 2.0;
			f[j] = fun(x[j]);
		}
		quanc8nofun = quanc8nofun + 8;
		double step = (x[16] - x0) / 16.0;
		double qleft = (w0*(f0 + f[8]) + w1*(f[1] + f[7]) + w2*(f[2] + f[6])
			+ w3*(f[3] + f[5]) + w4*f[4]) * step;
		qright[lev + 1] = (w0*(f[8] + f[16]) + w1*(f[9] + f[15]) + w2*(f[10] + f[14])
			+ w3*(f[11] + f[13]) + w4*f[12]) * step;
		double qnow = qleft + qright[lev + 1];
		double qdiff = qnow - qprev;
		area = area + qdiff;

		//  ***   stage 4 *** interval convergence test

		double esterr = fabs(qdiff) / 1023.0;
		if (abserr >= relerr*fabs(area))
			tolerr = abserr;
		else
			tolerr = relerr*fabs(area);
		tolerr = tolerr*(step / stone);

		// multiple logic conditions for the convergence test
		//key = 1;
		if (lev < levmin) key = 1;
		else if (lev >= levmax)
			key = 2;
		else if (quanc8nofun > nofin)
			key = 3;
		else if (esterr <= tolerr)
			key = 4;
		else
			key = 1;

		switch (key) {
			// case 1 ********************************* (mark 50)
		case 1:
			//      ***   stage 5   ***   no convergence
			//      locate next interval.
			nim = 2 * nim;
			lev = lev + 1;

			//      store right hand elements for future use.
			for (i = 1; i <= 8; i = i + 1)
			{
				fsave[i][lev] = f[i + 8];
				xsave[i][lev] = x[i + 8];
			}

			//      assemble left hand elements for immediate use.
			qprev = qleft;
			for (i = 1; i <= 8; i = i + 1)
			{
				j = -i;
				f[2 * j + 18] = f[j + 9];
				x[2 * j + 18] = x[j + 9];
			}
			continue;  // go to start of stage 3 "central calculation"
			break;

			// case 2 ********************************* (mark 62)
		case 2:
			quanc8flag = quanc8flag + 1.0;
			break;
			// case 3 ********************************* (mark 60)
		case 3:
			//    ***   stage 6   ***   trouble section
			//    number of function values is about to exceed limit.
			nofin = 2 * nofin;
			levmax = levout;
			quanc8flag = quanc8flag + (b - x0) / (b - a);
			break;
			// case 4 ********************************* (continue mark 70)
		case 4:
			break;
			// default ******************************** (continue mark 70)
		default:
			break;
			// end case section ***********************
		}

		//   ***   stage 7   ***   interval converged
		//   add contributions into running sums.

		quanc8result = quanc8result + qnow;
		quanc8errest = quanc8errest + esterr;
		cor11 = cor11 + qdiff / 1023.0;

		//  locate next interval

		while (nim != 2 * (nim / 2))
		{
			nim = nim / 2;
			lev = lev - 1;
		}
		nim = nim + 1;
		if (lev <= 0) break;  // may exit futher calculation

		//  assemble elements required for the next interval.

		qprev = qright[lev];
		x0 = x[16];
		f0 = f[16];
		for (i = 1; i <= 8; i = i + 1)
		{
			f[2 * i] = fsave[i][lev];
			x[2 * i] = xsave[i][lev];
		}
	}
	//  *** end stage 3 ***   central calculation

	//  ***   stage 8   ***   finalize and return

	quanc8result = quanc8result + cor11;

	//  make sure quanc8errest not less than roundoff level.

	if (quanc8errest == 0.0) return;
	do
	{
		temp = fabs(quanc8result) + quanc8errest;
		quanc8errest = 2.0*quanc8errest;
	} while (temp == fabs(quanc8result));

	return;
}

//# corresponding python code
//def sumg_c(z,u):
    //zplusu=z+u
    //zminusu=z-u
    //log_result = cmath.log((zplusu+1.0)/(zplusu-1.0))
    //plusresult = (1.0 -zplusu*zplusu)*log_result
    //log_result = cmath.log((zminusu+1.0)/(zminusu-1.0))
    //minusresult = (1.0 -zminusu*zminusu)*log_result
    //return(plusresult+minusresult)
    
inline dcomp sumg (dcomp z, dcomp u)
{
    dcomp zplusu=z+u;
    dcomp zminusu=z-u;
    dcomp log_result = log((zplusu + 1.0)/(zplusu - 1.0));
    dcomp plusresult = (1.0 - zplusu * zplusu) * log_result;
    log_result = log((zminusu + 1.0)/(zminusu - 1.0));
    dcomp minusresult = (1.0 - zminusu*zminusu)*log_result;
    return(plusresult+minusresult);
}


//inline dcomp g_c_quad (dcomp A)
//{   __complex128 QA;
    //__real__ QA = A.real();
    //__imag__ QA = A.imag();
    //__complex128 log_result = clogq((QA + 1.0)/(QA - 1.0));
    //__complex128 qout = (1.0 -QA*QA)*(log_result);
    //dcomp out= (dcomp) qout;
    //return out; 
//}


inline dcomp sumg_quad (dcomp z, dcomp u)
{    
    __complex128 zplusu, zminusu, qz,qu;
    __real__ qz = z.real();
    __imag__ qz = z.imag();
    __real__ qu = u.real();
    __imag__ qu = u.imag();
    
     zplusu=qz+qu;
     zminusu=qz-qu;
     __complex128 log_result = clogq((zplusu + 1.0)/(zplusu - 1.0));
     __complex128 plusresult = (1.0 - zplusu * zplusu) * log_result;
    log_result = clogq((zminusu + 1.0)/(zminusu - 1.0));
     __complex128 minusresult = (1.0 - zminusu*zminusu)*log_result;
     dcomp out=(dcomp)(plusresult+minusresult);
    return(out);
}
//#corresponding python code
//def Lindhard_LL(q,omega,gamma,w_p,U):  
    //w_c=complex(omega,gamma)
    //ww_used=  cmath.sqrt(w_c**2-complex(U,0.0)**2)
    //v_f=pow(w_p*w_p*3.0/4.0*math.pi,1.0/3.0)
    //prefactor=3*w_p*w_p/(v_f*v_f*q*q)
    //u= ww_used / (sqrt_2Q(q)*v_f)
    //z = sqrt_2Q(q) / (2 * v_f)
    //# d1 = g_c(z + u)
    //# d2 = g_c(z - u)
    //# f=0.5+ (d1 + d2)/(8*z) #following notation Sigmund's appendix and Lindhard 3.6
    //f=0.5+ sumg_c(z,u)/(8*z)
    //eps=1.0+prefactor*f
    //return eps


dcomp  Lindhard(double q, dcomp omega_c, double  omega0)
{   dcomp u, sumdterms;
    double z, Qrecoil;
   
    double v_f=pow(omega0*omega0*3.0/4.0*pi,1.0/3.0);
    
  
    if(Dispersion_relativistic && (q > 1.0))
    {   
        Qrecoil = sqrt(C * C * q * q + C * C * C * C) - C * C;
        z = sqrt(Qrecoil)/ (sqrt(2.0) * v_f); 
        u=omega_c/ (sqrt(2.0* Qrecoil) *v_f);
    }
    else
    {    Qrecoil=q*q/2.0;
         u=omega_c/ (q*v_f);
         z = q / (2 * v_f);
    }
	//double CHI_square = 1.0 / (pi*v_f);
    if (abs(u) < 500.0*z) // the transition momentum  may need some fine tuning
    {
       sumdterms=sumg(z,u);
    }
    else    
    {
        sumdterms=sumg_quad(z,u);
    }
    //double prefactor=3*omega0*omega0/(v_f*v_f*q*q);// wrong relativistically
  //  double prefactor=3*omega0*omega0/(4*z*z*pow(v_f,4));// this one works
    double prefactor=3*omega0*omega0/(2*Qrecoil*v_f*v_f);// this one works
    dcomp f = 0.5+sumdterms/(8*z);// note d1-d2 differs from Sigmund's f1+if2 (5.155) by 4z
    return (1.0+prefactor*f);
}

dcomp eps_LLX(double q, dcomp  omega_c, double omega0, double U)
{
    dcomp eps;
    dcomp omega_minus=sqrt(omega_c*omega_c-U*U);
    eps = Lindhard(q,omega_minus,omega0);
    return eps;
}

dcomp Lindhard_LL(double q, double  omega, double gamma, double  omega0, double omega_gap)
{
	dcomp z1, z2, z3, top, bottom, omega_c;
    omega_c=dcomp(omega, gamma);
    dcomp omega_minus=sqrt(omega_c*omega_c-omega_gap*omega_gap);
    if (q > q_transition*omega+0.00001)
    {
        double g_over_w = gamma / omega;
        z1 = dcomp(1.0, g_over_w);// omega should be unequal 0
        z2 =Lindhard(q,omega_minus,omega0) - dcomp(1.0, 0.0);
        if ( Apply_Mermin_Correction)
        {   
            dcomp omega_minus=dcomp(0.0,omega_gap+1e-10);   //else strange things happen when U=0
            z3 = Lindhard(q, omega_minus, omega0) - dcomp(1.0, 0.0);
            top = z1*z2;
            bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
            z1 = dcomp(1.0, 0.0) + top / bottom;
        }
        else
        {
            z1=  dcomp(1.0, 0.0) + z2;
        }
    }
    else //calculate the equivalent Drude-Lindhard
    {  
        if ( Apply_Mermin_Correction)
		{
			z1 = oneoverDL(q, omega, gamma, omega0, 1.0,omega_gap ); 
		}
		else  //plain lindhard, has double the nominal  width
		{
			z1 = oneoverDL(q, omega, 2* gamma, omega0,1.0,omega_gap );
		}
		
        z1 = dcomp(1.0, 0.0) / z1 ;
    }
    
	return (z1);
}

dcomp Lindhard_Direct(double q, double  omega, double gamma, double  omega0, double omega_gap)
{    dcomp eps;
     if (q > q_transition*omega+0.00001)
    {
        dcomp omega_c=dcomp(omega, gamma/2.0);
        dcomp ww_square_DL= omega_c*omega_c+ gamma*gamma/4.0-omega_gap*omega_gap;
        dcomp ww_used=sqrt( ww_square_DL);
        eps=Lindhard(q, ww_used,omega0);
    }
    else
    {
        dcomp z1 = oneoverDL(q, omega, gamma, omega0, 1.0,omega_gap );
        eps=dcomp(1.0, 0.0) / z1 ;
    }
    
	return eps;
}

dcomp   eps_Kaneko(double q, dcomp omega_c, double q_mean, double U, double gamma_fudge)
{   
    double gammaChi2=2*gamma_fudge/(q_mean);   //factor of 2 due to spin, original kaneko means gamma_fudge=0.5
    dcomp eps,dd,u;
    double z;
    dcomp omega_minus=sqrt(omega_c*omega_c-U*U);
   
    if(Dispersion_relativistic && (q > 1.0))
    {   
        double Qrecoil = sqrt(C * C * q * q + C * C * C * C) - C * C;
        u=omega_minus/ (sqrt(2.0* Qrecoil) * q_mean);
        z = sqrt(Qrecoil)/ (sqrt(2.0) * q_mean); 
    }  
    else
    {
         u=omega_minus/ (q*q_mean);
         z = q / (2 * q_mean);
    }
    double relerror=0.0;
    dcomp term1 = dcomp(0,1)/2.0* Faddeeva::w(u + z, relerror);
    dcomp term2 = dcomp(0,1)/2.0* Faddeeva::w(u - z, relerror);
    dd=1.0/(4*z)*(term2-term1);
    eps=1.0 + gammaChi2/(z*z)*(dd);
    return eps;
   
}
dcomp KanekoDirect (double q, double omega, double gamma, double  alpha, double U,  double gamma_fudge)
{
    dcomp omega_c=dcomp(omega, gamma/2.0);
	double Q =sqrt(1.0/alpha);
	dcomp  ww_square=omega_c*omega_c+omega_c.imag()*omega_c.imag()-U*U;  
    dcomp ww_used=sqrt(ww_square); 
    dcomp eps=eps_Kaneko(q, ww_used, Q, 0.0, gamma_fudge);
    return eps;
}



dcomp Kaneko(double q, double omega, double gamma, double  alpha, double U,  double gamma_fudge)
{
	dcomp z1, z2, z3, top, bottom, tmp;
	double g_over_w;
    double Q =sqrt(1.0/alpha);

    g_over_w = gamma / omega;
    dcomp omega_c= dcomp(omega,gamma);
    z1 = dcomp(1.0, g_over_w);// omega should be unequal 0 
    z2 = eps_Kaneko(q, omega_c, Q, U, gamma_fudge) ;
    
    if(!Apply_Mermin_Correction) return z2;
    z2 = z2 - dcomp(1.0, 0.0);
    z3 = eps_Kaneko(q, dcomp(0.0, 0.0), Q,U,  gamma_fudge) - dcomp(1.0, 0.0);
    top = z1*z2;
    bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
    z1 = dcomp(1.0, 0.0) + top / bottom;   
	return (z1);
}

dcomp diff_eps(double q, double omega, double gamma, double  alpha, double U, double  gamma_fudge)
{  
    dcomp eps_minus,eps_plus;
    double delta_alpha=0.005*alpha;
    if(!DirectMethod)// here we differentiate  after mermin correction if Apply_Mermin_Correction=true. not 100% sure about that one, seems to work
    {
		eps_minus = Kaneko(q, omega, gamma, alpha-delta_alpha,U,gamma_fudge);
		eps_plus  = Kaneko(q, omega, gamma, alpha+delta_alpha,U, gamma_fudge);
	}
	else
	{ 
		eps_minus = KanekoDirect(q, omega, gamma,  alpha-delta_alpha, U, gamma_fudge);
		eps_plus =  KanekoDirect(q, omega, gamma,  alpha+delta_alpha, U, gamma_fudge);
	}
    return (eps_plus-eps_minus)/(2.0*delta_alpha);
}
dcomp diff2_eps(double q, double omega, double gamma, double  alpha,double U,  double  gamma_fudge)
{  
    dcomp eps_minus,eps_plus;
    double delta_alpha=0.005*alpha;
    eps_minus = diff_eps(q, omega, gamma, alpha-delta_alpha,U,gamma_fudge);
    eps_plus  = diff_eps(q, omega, gamma, alpha+delta_alpha,U,gamma_fudge);
    return (eps_plus-eps_minus)/(2.0*delta_alpha);
}
dcomp diff3_eps(double q, double omega, double gamma, double  alpha, double U, double  gamma_fudge)
{  
    dcomp eps_minus,eps_plus;
    double delta_alpha=0.005*alpha;
    eps_minus = diff2_eps(q, omega, gamma,  alpha-delta_alpha, U, gamma_fudge);
    eps_plus  = diff2_eps(q, omega, gamma, alpha+delta_alpha, U, gamma_fudge);
    return (eps_plus-eps_minus)/(2.0*delta_alpha);
}





dcomp  calculate_chi_AA_LL(double q)  // Archubi-Arista suggested use of  Levine Louie to calculate Kaneko with gap
                                        // revert to plane Kaneko if energy edge  (=U)is 0, calculated at w_global
{
    dcomp eps,sum_chi, oneoverchi,sum_oneoverchi,  omega_minus, omega_c,z1;
    double alpha, Volume, Volumefraction,gamma_fudge;
    double ratio,U;

  
 
    sum_chi=dcomp(0.0,0.0);
    sum_oneoverchi=dcomp(0.0,0.0);
    for (int i = 0; i < MAXKANEKO; i++)
	{ 
		if ((fabs(N_Kaneko[i]) > 1E-10)) // && (w_global> Edge_ArchubiKaneko[i]) )
        {  
            U =  Edge_ArchubiKaneko[i];
            omega_c = dcomp( w_global, Width_Kaneko[i]);         
            alpha=1.0/(Q_Kaneko[i]*Q_Kaneko[i]);
            gamma_fudge=gamma_Kaneko[i];
         
          
            ratio=w_pl_l[i]*w_pl_l[i]* A_l_Kaneko[i]/(A_0_Kaneko[i]*w_pl_0[i]*w_pl_0[i]) ;    
            
            Volume=4.0*pi*N_Kaneko[i]/(w_pl_l[i]*w_pl_l[i]); 
            Volumefraction=Volume* UnitCellDensity; //gamma fudge already incorportaed via w_pl[i]
            if(Volumefraction > 1.0)
            {
                printf("warning, unphysical density Kaneko component  %i volumefraction: %6.4f\n",i, Volumefraction);
            }
            
   

            if (q < q_transition*w_global+0.00001)
            {
               //  printf("was here q= %6.4f\n",q);
                if ( Apply_Mermin_Correction||DirectMethod )
                {
                    z1 = oneoverDL(q, w_global, Width_Kaneko[i], w_pl_l[i], 1.0,U ); 
                }
                else  //plain RPA , has double the nominal  width
                { 
                    z1 = oneoverDL(q, w_global, 2* Width_Kaneko[i], w_pl_l[i], 1.0,U );
                }
                
                eps = dcomp(1.0, 0.0) / z1 ;
            }
          
            else if(l_Kaneko[i]==0)
            {   if(!DirectMethod) 
				{                     
					eps = Kaneko(q, w_global, Width_Kaneko[i], alpha, U, gamma_fudge);
			    }
			    else
			    { 
					 eps = KanekoDirect(q, w_global, Width_Kaneko[i],   alpha, U, gamma_fudge);
				}
            }
            else if  (l_Kaneko[i]==1)
            {  
                eps=  dcomp(1.0,0.0)-ratio*diff_eps(q, w_global, Width_Kaneko[i],   alpha,U, gamma_fudge);
            }
            else if  (l_Kaneko[i]==2)
            {   
                eps=  dcomp(1.0,0.0)+ratio*diff2_eps(q, w_global, Width_Kaneko[i],   alpha, U,  gamma_fudge);
            }
             else if  (l_Kaneko[i]==3)
            {   
                eps=   dcomp(1.0,0)-ratio*diff3_eps(q, w_global, Width_Kaneko[i],   alpha,U,  gamma_fudge);
            }
            dcomp oneoverchi=dcomp(1.0,0)/eps-dcomp(1.0,0.0);
            if(AddELF)
            {
                sum_oneoverchi+=oneoverchi*Volumefraction;
            }
            else
            {   eps=dcomp(1.0,0.0)/(Volumefraction*oneoverchi+dcomp(1.0,0.0));
                sum_chi+=(eps-dcomp(1.0,0.0));  
            }
        }
    }
    if(AddELF)
    {
        sum_chi=dcomp(1.0,0.0)/(sum_oneoverchi+dcomp(1.0,0.0))-dcomp(1.0,0.0);
    }
    return sum_chi;
}

double  calculate_Loss_AA_LL(double q)  // Archubi-Arista suggested use of  Levine Louie to calculate Kaneko with gap
                                        // revert to plane Kaneko if energy edge  (=U)is 0, calculated at w_global
                                    
{
    dcomp eps=calculate_chi_AA_LL(q)+dcomp(1.0,0.0);
    dcomp oneovereps= dcomp(1.0,0.0)/eps;
    double loss= -oneovereps.imag();
    return loss;
}                                        

dcomp  ChiDrude(double q, double current_w,  double  gamma, double w, double alpha,double  U)
{   dcomp chi, gamma_c; 
    double w_at_q_square, Qrecoil;
    gamma_c=dcomp(0,gamma+current_w/1000.0);   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    if(Dispersion_relativistic)
    {
        Qrecoil = sqrt(C * C * q * q + C * C * C * C) - C * C;  //m_e is 1 in a.u.
    }
    else
    {
        Qrecoil = q * q / 2.0;
    }   
    if (Full_dispersion)
    {
        double k_f=pow(w*w*3.0/4.0*pi,1.0/3.0);
        w_at_q_square = w*w + 2.0 * alpha *k_f*k_f/3.0 *Qrecoil +Qrecoil * Qrecoil;// Lundquist way, alpha should be 1 2 because Q=q^2/2
    }
    else
    {
        // now there are different alpha's in use in the literature here we use dispersion relative to free particle
        w_at_q_square = (w + alpha * Qrecoil)*(w + alpha * Qrecoil);
    }
    chi =  -1.0 / (w_global*w_global - w_at_q_square - U*U + w_global*gamma_c);// U tries to mimick band gap here, just as in LL
    return chi;
}



dcomp DrudeSum(double q)
{
	dcomp  oneovereps,chi_this_osc, eps, one_over_chi_total;
	one_over_chi_total=dcomp(0.0,0.0);
	eps=dcomp( epsbkg,0.0);
    oneovereps=dcomp(1.0/epsbkg,0.0);

	for (int i = 0; i < MAXOSC; i++)
	{
		if (fabs(Ai[i]) >1e-90)
		{  
            chi_this_osc=Ai[i]*ChiDrude(q, w_global, gammai[i], wi[i], alphai[i],gapi[i]);
            if(AddELF)
            {
			    oneovereps+=(dcomp(1.0,0.0)/(dcomp(1.0,0.0)+chi_this_osc))-dcomp(1.0,0.0);
            }
            else
            {
                eps+=chi_this_osc;
            }
		}
	}
    if(AddELF)
    {  
        eps= (dcomp(1.0,0.0)/oneovereps);
    }
 
    return (eps);
}

dcomp  oneoverDL(double q, double current_w, double  gamma, double w, double alpha, double U)
{
	double w_at_q_square, Qrecoil;

    if(Dispersion_relativistic)
    {
        Qrecoil = sqrt(C * C * q * q + C * C * C * C) - C * C;  //m_e is 1 in a.u.
    }
    else
    {
        Qrecoil = q * q / 2.0;
    }     
                  
    if (Full_dispersion)
    {  
        double k_f=pow(w*w*3.0/4.0*pi,1.0/3.0);
        w_at_q_square = w*w + 2.0 * alpha *k_f*k_f /3.0 *Qrecoil +Qrecoil * Qrecoil;  //Lundqvist way, alpha should be 1, 
    }
    else
    {
        // now there are different aplha's in use in the literature here we use dispersion relative to free particle
        w_at_q_square = pow(w + alpha * Qrecoil,2);
    }
    dcomp gamma_C(0.0,gamma+current_w/1000.0);   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    dcomp one_over_eps = 1.0 + w*w/(current_w*current_w - w_at_q_square - U*U+ current_w*gamma_C );
	return one_over_eps;
}
double TaucNormalisation(double  C, double E0_thisq,  double Egap)
{   int nstep=1000;
    double stepsize= C/40;
    double sum=0.0; 
    double x;
    double Elower= E0_thisq;
    double Eupper= E0_thisq;
    
    for (int i=0; i < nstep/2; i++)
    {
        if( Elower > Egap)
        {
            x = E0_thisq*C*pow((Elower-Egap),2)/((pow((Elower*Elower-E0_thisq*E0_thisq),2)+C*C*Elower*Elower))*stepsize;
            sum +=x;
        }
        x = E0_thisq*C*pow((Eupper-Egap),2)/((pow((Eupper*Eupper-E0_thisq*E0_thisq),2)+C*C*Eupper*Eupper))*stepsize;
        sum +=x;
        stepsize=stepsize*(1.07);
        Eupper += stepsize;
        Elower -= stepsize;
        //if(Eupper > (E0_thisq + 50*(C +5*C*C/(E0_thisq-Egap)))) //if E0 close to gap peak becomes wide
        //{  // printf("normalisation ended n=%i\n",i);
            //return(sum);
        //}
        if(x/sum < 1e-6) 
        {   printf("last i %i, sum %6.4f\n",i,sum);
            return (sum);
        }
            
    } 
   // printf("sum  %6.4f\n",sum);
    return sum;
}
dcomp  ChiTaucLorentz(double E, double  C, double E0,  double Egap)
{
    double Chi1,Chi2;
    double a_ln=(Egap*Egap-E0*E0)*E*E+Egap*Egap*C*C-E0*E0*(E0*E0+3*Egap*Egap);
    double a_atan=(E*E-E0*E0)*(E0*E0+Egap*Egap)+Egap*Egap*C*C;
    double gamma=sqrt(E0*E0-C*C/2);
    double alpha=sqrt(4*E0*E0-C*C);
    double zeta4=(E*E-gamma*gamma)*(E*E-gamma*gamma)+alpha*alpha*C*C/4;
    
    if(E> Egap)
    {
        Chi2=E0*C*pow((E-Egap),2)/((pow((E*E-E0*E0),2)+C*C*E*E)*E);
       
    }
    else
    {   
        Chi2=0.0;
    }
        
    Chi1  = C/(pi*zeta4)*a_ln/(2*alpha*E0)*log((E0*E0+Egap*Egap+alpha*Egap)/(E0*E0+Egap*Egap-alpha*Egap));
    Chi1 -= a_atan/(pi*zeta4*E0)*(pi-atan((2*Egap+alpha)/C)+ atan((-2*Egap+alpha)/C));
    Chi1 += 2*E0/(pi*zeta4*alpha)*Egap*(E*E-gamma*gamma)*(pi+2*atan(2*(gamma*gamma-Egap*Egap)/(alpha*C)));
    Chi1 -= E0*C/(pi*zeta4)*(E*E+Egap*Egap)/E*log(fabs(E-Egap)/(E+Egap));
    Chi1 += 2*E0*C/(pi*zeta4)*Egap*log((fabs(E-Egap)*(E+Egap))/sqrt(pow((E0*E0-Egap*Egap),2)+Egap*Egap*C*C));

    return dcomp(Chi1,Chi2);
}


dcomp TaucLorentz_Sum(double q)
{          
    dcomp ChiTL, sumoneovereps,sumeps,sumoneoverchi;
    // first the shift due to dispersion
    double Qrecoil,currentE0;
    double norm_this_q;
    if(Dispersion_relativistic)
    {
        Qrecoil = sqrt(C * C * q * q + C * C * C * C) - C * C;  //m_e is 1 in a.u.
    }
    else
    {
        Qrecoil = q * q / 2.0;
    }     
         
	sumoneoverchi = dcomp(0.0, 0.0);
    sumeps=dcomp(epsbkg,0.0);   
    
	for (int i = 0; i < MAXOSC; i++)
	{
		if (fabs(Ai[i]) >1e-90)
		{   currentE0=wi[i]+alphai[i]*Qrecoil;
            norm_this_q=TaucNormalisation(gammai[i], currentE0,  gapi[i]);
            ChiTL = (Ai[i]*TaucNormalisationq0[i]/norm_this_q)*ChiTaucLorentz(w_global,gammai[i],currentE0, gapi[i]);
            sumeps+=ChiTL ;
        }
    }
    return (sumeps);
}

dcomp DL_Sum(double q)
{
	dcomp oneoverepsDL, sumoneovereps,sumeps,epsDL,sumoneoverchi;
	sumoneoverchi = dcomp(0.0, 0.0);
    sumeps=dcomp(1.0,0.0);
    
	for (int i = 0; i < MAXOSC; i++)
	{
		if (fabs(Ai[i]) >1e-90)
		{       
            oneoverepsDL = oneoverDL(q,w_global,gammai[i], wi[i], alphai[i], gapi[i]);
            if(AddELF)
            {
                sumoneoverchi+= + Ai[i] * (oneoverepsDL-dcomp(1.0,0.0));  //add "1/Chi"'s not 1/eps
            }
            else
            {   epsDL=dcomp(1.0,0.0)/(Ai[i]*(oneoverepsDL-dcomp(1.0,0.0))+dcomp(1.0,0.0));
                sumeps+=(epsDL-dcomp(1.0,0.0));
                
            }
		}
	}
    if(AddELF)
    {
	    sumoneovereps=sumoneoverchi +dcomp(1.0,0.0);
	    return(dcomp(1.0, 0.0) / sumoneovereps);
    }
    else
    {
        return (sumeps);
    }
}

dcomp LindhardLL_Sum(double q)
	{   
    dcomp epsMerminLL,oneoverepsMerminLL, sumoneovereps,sumeps,sumoneoverchi, epsMLL;
    double gamma_used;
	sumoneoverchi = dcomp(0.0, 0.0);
    sumeps=dcomp(1.0,0.0);

	for (int i = 0; i < MAXOSC; i++)
	{
		if (fabs(Ai[i]) > 0.0)
		{    
			if (abs((w_global - wi[i]) / gammai[i]) < 100000)  
			{   
				gamma_used=gammai[i];
			}
			else
			{ 
				gamma_used=(w_global - wi[i])/10000;
			}
			if(!DirectMethod)
			{
				epsMerminLL = Lindhard_LL(q, w_global, gamma_used, wi[i], gapi[i]);
			}
			else  //Direct obviously	
			{
				epsMerminLL = Lindhard_Direct(q, w_global, gamma_used, wi[i], gapi[i]);
			}
            oneoverepsMerminLL = dcomp(1.0,0.0)/epsMerminLL ;
            if(AddELF)
            {
                sumoneoverchi+= + Ai[i] * (oneoverepsMerminLL-dcomp(1.0,0.0));  //add "1/Chi"'s not 1/eps
            }
            else
            {   
                epsMLL=dcomp(1.0,0.0)/(Ai[i]*(oneoverepsMerminLL-dcomp(1.0,0.0))+dcomp(1.0,0.0));
                sumeps+=(epsMLL-dcomp(1.0,0.0));
                
            }
		}
	}
    if(AddELF)
    {
	    sumoneovereps=sumoneoverchi +dcomp(1.0,0.0);
	    return(dcomp(1.0, 0.0) / sumoneovereps);
    }
    else
    {
        return (sumeps);
    }
}

dcomp calculate_eps_Belkacem(double q, double omega,double A_i, double w_i, double gamma_i)

{     
    dcomp eps, chi, gamma, fraction1,fraction2;
    double W_plasmon_Belkacem;
    double omega_k,C1,M1,currentPoisson, ratio;//, w_i_dens;
    int n,n_min, n_max;
    
     if(q> q_transition)
    {
        if(Dispersion_relativistic)
        {
            omega_k = sqrt(C * C * q * q + C * C * C * C) - C * C;  //m_e is 1 in a.u.
        }
        else
        {
            omega_k = q * q / 2.0;
        }
    }
    else  // use a minimum value for q
          //for small q all intensity in n=1 component and then omega_k cancels
    {
        omega_k=q_transition*q_transition/2.0;
    }
    ratio=omega_k/w_i;
    W_plasmon_Belkacem=PL_GOS*sqrt(A_i); //PL_GOS plasmon energy of 1 electron per unit cell. Ai_Belkacem is number of oscillators per unit cell
  //  w_i_dens = sqrt(w_i*w_i+W_plasmon_Belkacem*W_plasmon_Belkacem);//density effect may affect os. frequency  
    n_min=round(omega/w_i)-16;  // was 60
    if (n_min < 1) n_min=1;
    n_max=round(omega/w_i) + 16;  // was 60 divided w_i_dens
    double fraction_covered= Poisson( 0, ratio);
    chi=dcomp(0.0,0.0);
    gamma = dcomp(0.0, gamma_i);//+ omega/10000);
    C1=pow(W_plasmon_Belkacem,2)/(2.0*omega_k);
    for(n=n_min; n< n_max; n++)
    {    
         currentPoisson= Poisson( n, ratio);
         fraction_covered += currentPoisson;
         M1=C1* currentPoisson;
         fraction1= 1.0/(n*w_i - omega - gamma);
         fraction2= 1.0/(n*w_i + omega + gamma);
         chi += M1*(fraction1 + fraction2);
    }
  
   
    
    eps=dcomp(1.0,0)+chi;
    return eps;
}

 
dcomp calculate_chi_Belkacem(double q)  //calculate at w_global
{   
    dcomp  chi, sum_chi,sum_one_over_eps,one_over_eps;
    sum_chi=dcomp(0.0,0.0);
    sum_one_over_eps=dcomp(1.0,0.0);
    //in version 10 we applied Merminization to Belkacem. removed, does not make sense to me anymore
    for (int i = 0; i < MAXBELKACEM; i++)
    { 
        if (Ai_Belkacem[i] > 1E-10)
        {
             chi=  calculate_eps_Belkacem(q, w_global,Ai_Belkacem[i], wi_Belkacem[i],gammai_Belkacem[i])-dcomp(1.0,0);
             if(AddELF)
             {
                one_over_eps=dcomp(1.0,0.0)/(dcomp(1.0,0.0)+chi);
                sum_one_over_eps+=one_over_eps-dcomp(1.0,0.0);  //add  one over chi's
             }
             else
             {
                 sum_chi+=chi;
             }
        }    
    }
    if(AddELF)
    { 
        sum_chi= dcomp(1.0,0.0)/sum_one_over_eps-dcomp(1.0,0);
    }
    return(sum_chi);

 } 
   
double calculate_loss_Belkacem(double q)  //calculate at w_global
{   dcomp chi_Belk,oneovereps;      
    chi_Belk=calculate_chi_Belkacem( q);
    oneovereps=dcomp(1.0,0.0)/(chi_Belk+dcomp(1.0,0.0));
    return -oneovereps.imag();
}

dcomp calculate_eps_osc(double q)
{
	dcomp eps;
    eps=dcomp(1.0,0.0);// so initialized if no oscillators
	if ((modelMerminLL) )
	{
		eps = LindhardLL_Sum(q);
	}
	else if (modelDrude)
	{   
		eps = DrudeSum(q);
	}
	else if (modelDL)
	{   
		eps = DL_Sum(q);
	}
    else if (modelTaucLorentz)
	{   
		eps = TaucLorentz_Sum(q);
	}
   	return eps;
}

dcomp get_eps_component(int j,double q)
{   
    dcomp eps, gamma;
    gamma = dcomp(0.0, gammai[j]);
  
    if (modelMerminLL)
    {
        eps = Lindhard_LL(q, w_global, gammai[j], wi[j], gapi[j]);
    }
    else if (modelDrude)
    {
        eps = dcomp(1.0,0)+  Ai[j]*ChiDrude(q,w_global,gammai[j], wi[j], alphai[j], gapi[j]);
    }
    else if (modelDL)
    {   
        eps =dcomp (1.0,0)/oneoverDL(q,w_global,gammai[j], wi[j], alphai[j], gapi[j]);
    }
  
    return eps;
}


dcomp calculate_chi_GOS_dens(double q, double omega, int GOSi) //// calculates eps1 eps2 from GOS including density effect
{  
    double chi_real,chi_imag;
    double W_plasmon_GOS_square;
    double  currentscaling;
 
	int iq  = int(q/GOSScalingStepSize);  // 
	//double q_remainder = fmod(q, GOSScalingStepSize);
    double position_in_bin=(q-iq*GOSScalingStepSize)/GOSScalingStepSize;
    if (ApplySumRuleToGOS)
			{   double before,after;
				if (iq > 79)  iq = 79;
				before = GOS_ScalingFactor[GOSi][iq];
				if (iq < 79) after = GOS_ScalingFactor[GOSi][iq + 1];
				else after = before;
				currentscaling = before*(1.0 - position_in_bin) + position_in_bin*after;//currentscaling so we obey sum rule
			
			}
			else 
			{   
				currentscaling=1.0;
			}

    W_plasmon_GOS_square=PL_GOS*PL_GOS*AiGOS[GOSi]*currentscaling; //PL_GOS plasmon energy of 1 electron per unit cell. Ai is number of GOS electrons per unit cell of comp. GOSi
                                                                   // contains the scaling factor
    chi_real=0.0;
    chi_imag=pi*W_plasmon_GOS_square/(2.0*omega)*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, omega);
   
    if(omega < maxEnergyDensityEffect)
    {
      
        double current_StepSize=Edgei_GOS[GOSi]/200;// assume these are core levels  energy resolution required of the order of 0.5% of BE?
        if(current_StepSize < StepSize) current_StepSize= StepSize;
        double w_below=omega-current_StepSize;
        double w_above=omega+current_StepSize;
        int KKWidth=5000;
        for(int i=1; i< KKWidth; i++)
        {
            w_below-=current_StepSize;
            w_above+=current_StepSize;
            double Contribution_top=0.0;
            double Contribution_below=0.0;
            if(w_below >StepSize/2)
            {   
                Contribution_below=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_below)/(omega*omega- w_below*w_below);
            }
            Contribution_top=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_above)/(omega*omega- w_above*w_above);
            chi_real-=(Contribution_below+Contribution_top);
         
            if (kk_variable_stepsize)current_StepSize=current_StepSize*1.01;
        }
        
    }
    else //too small to bother
    { 
        chi_real=0.0;
    }
    return dcomp(chi_real,chi_imag);
}

double calculate_loss_GOS_dens(double q)
{ // calculate contribution GOS to loss function including density effect
    dcomp oneovereps, sum_chi, chiGOS;
    double loss=0.0;
    sum_chi=dcomp(0.0,0.0);
    	for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0.0)
            {
                chiGOS=calculate_chi_GOS_dens(q, w_global,i);
                if(AddELF)
                {
                    oneovereps=dcomp(1.0,0.0)/(chiGOS+dcomp(1.0,0.0));
                    loss-=oneovereps.imag();
                }
                else
                {
                    sum_chi+=chiGOS;
                }
                
            }
        }
        if(!AddELF)
        {
            oneovereps=dcomp(1.0,0.0)/(sum_chi+dcomp(1.0,0.0));
            loss=-oneovereps.imag();
        }
     return loss;  
}
dcomp calculate_chi_allGOS_dens( double q)
{ // calculate eps corresponding to the GOS, incl density effect
    dcomp chiGOS,chiGOS_i, oneoverchi, sum_chi, sum_oneoverchi;
    sum_chi=dcomp(0.0,0.0);
    sum_oneoverchi=dcomp(0.0,0.0);
    	for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0.0)
            {  if (errno !=0) my_perror("before chi_allGOS an error occured");
              chiGOS_i=calculate_chi_GOS_dens(q, w_global,i);
             
                if(AddELF)
                {
                    oneoverchi=dcomp(1.0,0.0)/(chiGOS_i+dcomp(1.0,0.0))-dcomp(1.0,0.0);
                    sum_oneoverchi+=oneoverchi;
                }
                else
                {
                    sum_chi+=chiGOS_i;
                }
            }
            
        }
        if(AddELF)
        {
            chiGOS=dcomp(1.0,0.0)/(sum_oneoverchi+dcomp(1.0,0.0)) -dcomp(1.0,0.0);
        }
        else
        {
            chiGOS=sum_chi;
        }
     return chiGOS;  
            
}



double  mylossfun(double q)
{
	dcomp eps;
	double result,GOS,Kaneko,Belkacem;
	eps = calculate_eps_osc( q);

	result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
    
	//evaluateGOS
    GOS=calculate_loss_GOS_dens(q);
  
   
    Kaneko= calculate_Loss_AA_LL(q);
    Belkacem=calculate_loss_Belkacem(q);
    
    result = result + GOS + Kaneko + Belkacem;
   
	return(result/q);  //return epsilon over q, as required for quanc8 integral
}



double  fun_neutral(double q)  // have to add Belkacem and Kaneko?
{
//	Just change Z ^ 2 -> (Z - rho(q)) ^ 2
//		for H0 it reads
//			Z = 1
//			rho(q) = 1 / (1.0 + (2 * q) ^ 2) ^ 4
	dcomp eps;
	double result, GOS,rho, factor;
	eps = calculate_eps_osc(q);

	result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares

    //GOS = calculate_GOS(q);
    GOS=calculate_loss_GOS_dens(q);
    result = result + GOS;
	rho = 1.0 / ((1.0 +  q*q /4.0)*(1.0 +  q*q/4.0));  //see email Pedro March 1 2019
	factor = (1.0 - rho)*(1.0 - rho);
	return(factor*result / q);  //return epsilon over q, as required for quanc8 integral
}



double DSEPfun(double q)
{  //for any incoming direction
	dcomp eps, epsfraction;
	double qs1_plus,qs1_minus, result;
   
    double theta_E; // characteristic angle Egerton EELS (3rd ed. ) eq. 3.28 omega/ gamma E_0
    double theta_scat,tmp;   //calculated here in the small angle approximation
    
    theta_E= (w_global/E_0)*  Egerton_rel_cor_factor;
    tmp=q*q-p_0*p_0*theta_E*theta_E;
    if(tmp < 0.0) 
    {
    // |tmp| values are very small (1e-14) so contrinute nothing
        return 0.0;
    }
    theta_scat=sqrt(tmp)/p_0;
	if (modelMerminLL)
	{
		eps = LindhardLL_Sum(q);
	}
	else if (modelDrude)
	{
		eps = DrudeSum(q);
	}
	else if (modelDL)
	{
		eps = DL_Sum(q);
	}
    
	epsfraction = (eps - dcomp(1.0, 0.0))*(eps - dcomp(1.0, 0.0)) / (eps*(eps + dcomp(1.0, 0.0)));
 
    qs1_plus=abs(p_0*theta_scat*cos(theta_global)+p_0* theta_E * sin (theta_global));
    qs1_minus=abs(p_0*theta_scat*cos(theta_global)-p_0 * theta_E * sin (theta_global));
    result = epsfraction.imag() *(qs1_plus+ qs1_minus )/ (q*q*q);
    return (result);
}

double DSEPFun_per_Os(double q)

{  //for any incoming direction
	dcomp eps, epsfraction;
	double qs1_plus,qs1_minus, result;
   
    double theta_E; // characteristic angle Egerton EELS (3rd ed. ) eq. 3.28 omega/ gamma E_0
    double theta_scat,tmp;   //calculated here in the small angle approximation
    
    theta_E= (w_global/E_0)*  Egerton_rel_cor_factor;
    tmp=q*q-p_0*p_0*theta_E*theta_E;
    if(tmp < 0.0) printf("problems\n");
    theta_scat=sqrt(tmp)/p_0;
    eps =  get_eps_component(I_Os_global,q);
    
	epsfraction = (eps - dcomp(1, 0))*(eps - dcomp(1, 0)) / (eps*(eps + dcomp(1, 0)));
    qs1_plus=abs(p_0*theta_scat*cos(theta_global)+p_0* theta_E * sin (theta_global));
    qs1_minus=abs(p_0*theta_scat*cos(theta_global)-p_0 * theta_E * sin (theta_global));
    result = epsfraction.imag() *(qs1_plus+ qs1_minus )/ (q*q*q);
    return (result);
}

     dcomp calc_total_eps( double omega, double q)
// the addition of the different contributions is an approximation only and not valid at significant densities of more than one type
// problems can be expected if there are  contributions of different models at similar energies
	{  
        dcomp eps,chi_Belk, chi_GOS,chi_Kaneko, totaleps, oneovereps; 
        dcomp oneover_chi_GOS,oneover_chi_Belk,oneover_chi_Kaneko;
        w_global = omega;
        eps=calculate_eps_osc(q);
        chi_Belk=calculate_chi_Belkacem( q);
        chi_GOS=calculate_chi_allGOS_dens(q);
        chi_Kaneko=calculate_chi_AA_LL(q);
        if (AddELF)
        {   
            oneovereps= dcomp(1.0,0.0)/eps; //oscillator contribution
            oneover_chi_Belk   =  dcomp(1.0,0.0)/(chi_Belk   + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            oneover_chi_GOS    =  dcomp(1.0,0.0)/(chi_GOS    + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            oneover_chi_Kaneko =  dcomp(1.0,0.0)/(chi_Kaneko + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            totaleps= dcomp(1.0,0.0)/(oneovereps+oneover_chi_Belk+oneover_chi_GOS+oneover_chi_Kaneko);
        }
        else
        {
            totaleps=eps+chi_GOS+chi_Belk+chi_Kaneko;
        }
        if (errno !=0) my_perror("calc_total_eps an error occured");
        return totaleps;
    }
    
    
    double Surface_loss_only(double omega, double q)
    {  
	dcomp eps, epsfraction;
    eps = calc_total_eps(omega, q);
	epsfraction = (eps - dcomp(1.0, 0.0))*(eps - dcomp(1.0, 0.0)) / (eps*(eps + dcomp(1.0, 0.0)));
	return (epsfraction.imag());  
}

int  copyP_to_Vars(double *p,  int modelchoice)
{

	double SumAi;
    double precision;
    int Param_Offset;
	modelMerminLL = false; modelDrude = false; modelDL = false; DirectMethod=false;modelTaucLorentz = false;
	if (modelchoice == 1)
	{
		 modelDrude = true;
	}
	else if (modelchoice == 2)
	{
		modelDL = true;
	}
	else if (modelchoice == 3)
	{
		modelMerminLL = true; 
	}
    else if (modelchoice == 4)
	{
		modelTaucLorentz = true; 
	}
    epsbkg=1.0;
	if (modelTaucLorentz) epsbkg = p[0]; 
   

	// now read the components
	SumAi = 0.0; 

	for (int i = 0; i < MAXOSC; i++)
	{  
		Ai[i] = p[5* i + 1];
		wi[i] = p[5 * i + 2] / Hartree;
		gammai[i] = p[5 * i + 3] / Hartree;
        if(fabs(Ai[i])> 1e-90 )
        {
            if (modelDrude) Ai[i] = Ai[i] / (Hartree*Hartree);// Ai is in eV^2 in DL
            alphai[i] = p[5 * i + 4];
            gapi[i] = p[5* i + 5] / Hartree;  //  gapi in LL model,gapi[]=0 gives Mermin
            if (modelTaucLorentz)
            {
                 Ai[i] = Ai[i] / (Hartree);// Ai is in eV in DT
                 TaucNormalisationq0[i]=TaucNormalisation(gammai[i],wi[i],gapi[i]);
            }  
            SumAi = SumAi + Ai[i];
        }
	}
    Param_Offset = 5* MAXOSC;
    for (int i = 0; i < MAXGOS; i++)
    {   
        AiGOS[i] = p[4 *  i+ Param_Offset + 1];
        Zi_GOS[i] = p[4 * i + Param_Offset + 2];
        int nl = (int)round(p[4 * i + Param_Offset + 3]);
        if (nl == 10)
        {
            n_i_GOS[i ] = 1;
            l_i_GOS[i ] = 0;
        }
        else if (nl == 20)
        {
            n_i_GOS[i] = 2;
            l_i_GOS[i] = 0;
        }
        else if (nl == 21)
        {
            n_i_GOS[i ] = 2;
            l_i_GOS[i ] = 1;
        }
        else if (nl == 30)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i ] = 0;
        }
        else if (nl == 31)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i] = 1;
        }
        else if (nl == 32)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i] = 2;
        }
        else
        {
            n_i_GOS[i] = 0;  //this should never occur
            l_i_GOS[i] = 0;
        }
        Edgei_GOS[i] = p[4 * i + Param_Offset + 4] / Hartree;

    }
    Param_Offset += 4 * MAXGOS;
     for (int i = 0; i < MAXBELKACEM; i++)
    {
        Ai_Belkacem[i]      = p[4 * i + Param_Offset + 1];
        wi_Belkacem[i]      = p[4 * i + Param_Offset + 2]/Hartree;
        gammai_Belkacem[i]  = p[4 * i + Param_Offset + 3]/Hartree;
    }
    Param_Offset += 4 *  MAXBELKACEM;
	

    for (int i = 0; i < MAXKANEKO; i++)
    {
        N_Kaneko[i]        = p[6 * i + Param_Offset + 1];
        Q_Kaneko[i]           = p[6 * i + Param_Offset + 2];
        Width_Kaneko[i]       = p[6 * i + Param_Offset + 3]/Hartree;
        Edge_ArchubiKaneko[i] = p[6 * i + Param_Offset + 4]/Hartree;
        l_Kaneko[i]= (int)round(p[6 * i + Param_Offset + 5]);
        gamma_Kaneko[i]       = p[6 * i + Param_Offset + 6];
        int l=l_Kaneko[i];
        //double fraction_of_shell_filled=N_Kaneko[i]/(2.0*(2*l+1));
        
        double Q=Q_Kaneko[i];
        double alpha = 1/(Q*Q);
        
        if((fabs(N_Kaneko[i]) > 1E-10)&&(!OriginalKaneko) )
            {  
            int  doublefact;
            if(l==0) doublefact=1;
            else if (l==1) doublefact=3;
            else if (l==2) doublefact=5*3;
            else  if (l==3) doublefact=7*5*3;
            else
            {
                printf(" lvalue %i not implemented for Kaneko \n",l);
                doublefact=0;  // to avoid warning used uninitialized
            }
           
            A_l_Kaneko[i]=pow(2.0,l)* pow(alpha,(2.0*l+3.0)/2.0)/(doublefact* pow(pi,1.5));
            A_0_Kaneko[i]=pow(alpha/pi,1.5);

           w_pl_l[i]=sqrt(gamma_Kaneko[i]*2*(2*l+1)*doublefact*pow(Q,3)*exp(l)/(pow(2.0*l,l)*sqrt(pi)) );// does not depend on occupation level, assumes all levels full
           w_pl_0[i]=sqrt(gamma_Kaneko[i]*2*pow(Q,3)/sqrt(pi));
        }
        else if((fabs(N_Kaneko[i]) > 1E-10)&&(OriginalKaneko) )
        {   Q_Kaneko[i]= Q_Kaneko[i]*pow(N_Kaneko[i],1.0/3.0);// variable Q_kaneko now contains Kaneko's  q_mean
            w_pl_0[i]=sqrt(pow(Q_Kaneko[i],3)/(sqrt(pi)) );  // no spin maximum occupation = 1
            w_pl_l[i]=w_pl_0[i]; // implicitly assumes l=0
            l_Kaneko[i]=0;   // implicitly assumes l=0
            gamma_Kaneko[i]=0.5* gamma_Kaneko[i];
         
           // Q_Kaneko[i]= Q_Kaneko[i]*pow(N_Kaneko[i],1.0/3.0);// variable Q_kaneko now contains Kaneko's  q_mean
           // printf("original Kaneko, plasmon energy %6.3f q_mean %6.2f\n",w_pl_l[i]*Hartree, Q_Kaneko[i] );
        }
    }

        

    Param_Offset += 6*MAXKANEKO;
	E_0= p[Param_Offset + 1]/Hartree;

	
	precision = p[Param_Offset + 2];
    abserr = 1e-50;
	//relerr = 1.0e-5;


    relerr = 2.0e-3 / precision;

	q2maxFactor = p[Param_Offset + 3];
	q2surfmaxFactor = p[Param_Offset + 4];
	NStep = (int)p[Param_Offset + 5];
	FirstEnergy = p[Param_Offset + 6] / Hartree;
	StepSize_eV = p[Param_Offset + 7];
    StepSize = StepSize_eV / Hartree;
    
       
    
    UnitCellDensity = p[Param_Offset + 8];// note this index is out of order, late addition  per ^3
	UnitCellDensity = UnitCellDensity*BohrRadius*BohrRadius*BohrRadius;// now per a.u.^3
    PL_GOS = sqrt(4.0*pi*UnitCellDensity);
    lin_cont_deltaE=p[Param_Offset + 9];
   
    q_lower=p[Param_Offset + 10];
    q_upper=p[Param_Offset + 11];
    LastMomentum=p[Param_Offset + 12];
    Stepsize_qplot=p[Param_Offset + 13];
    q_transition=p[Param_Offset + 14];  //where small q limit is reached (Mermin, Belkacem)
    ExchangeCorrection=p[Param_Offset + 15];
    if(p[Param_Offset + 16]==0.0)
    {
        OriginalKaneko=false;
    }
    else
    {
        OriginalKaneko=true;
    }
    if(p[Param_Offset + 17]==0.0)
    {
        AddELF=false;
    }
    else
    {
        AddELF=true;
    }
    if(p[Param_Offset + 18]==0.0)
    {
        ApplySumRuleToGOS=false;
    }
    else
    {
        ApplySumRuleToGOS=true;
    }
    
    if(p[Param_Offset + 19]==0.0)
    {
        Apply_Mermin_Correction=false;
        DirectMethod=false;
    }
    else if (p[Param_Offset + 19]==1.0)
    {
        Apply_Mermin_Correction=true;
        DirectMethod=false;
    }
    else
    {   Apply_Mermin_Correction=false;
        DirectMethod=true;
	}
    if(p[Param_Offset + 20]==0.0)
    {
        Full_dispersion=false;
    }
    else
    {
        Full_dispersion=true;
    }
    maxEnergyDensityEffect = p[Param_Offset + 21]/Hartree;
    if( maxEnergyDensityEffect < 0.0)
    {
        maxEnergyDensityEffect=-maxEnergyDensityEffect;
        kk_variable_stepsize=false;
        printf(" no variable stepsize\n");
    }
    else
    {
        kk_variable_stepsize=true;
        printf("using  variable stepsize\n");
    }    
   
    // for REELS spectrum 
    sigma=p[Param_Offset + 22]/2.355;
    PIcoef1= p[Param_Offset + 23];
    PIcoef2= p[Param_Offset + 24];
    PIcoef3= p[Param_Offset + 25];
    theta0= p[Param_Offset + 26]*pi/180.0;  //now in rad
    theta1= p[Param_Offset + 27]*pi/180.0;
    surf_ex_factor= p[Param_Offset + 28];
    if(p[Param_Offset + 29]==0.0)
    {
        Local_Desep=false;
    }
    else
    {
        Local_Desep=true;
    }
    fraction_DIIMFP= p[Param_Offset + 30];
    if(p[Param_Offset + 31]==1.0)
    { 
        proton=true;
        Projectile_Mass=MassProton;
        ExchangeCorrection=0;  // make sure  we are not bothered by exchange correction for protons
    }
    else
    {   proton=false;  //electron
        Projectile_Mass=1.0;
    }
    rest_mass_energy=Projectile_Mass*C*C;
    if(p[Param_Offset + 32]==1.0)
    {
        Dispersion_relativistic= true;
        printf("rel\n");
        
    }
    else
    {
         Dispersion_relativistic= false;
         printf("non_rel\n");
    }    
    theta_max= p[Param_Offset + 33]/1000.0;  // now in rad
    MottCorrection = false;
    if(p[Param_Offset + 34]==1.0 && proton) MottCorrection = true;
    
    Egerton_rel_cor_factor = (E_0+Projectile_Mass*C*C) / (E_0 +2.0*Projectile_Mass*C*C);// Egerton 3rd edition page 427 for theta_E
    rel_cor_factor =  (1 + E_0 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_0 / (Projectile_Mass*C*C)),2);	
    v_0 = sqrt(2.0*  rel_cor_factor*E_0/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    gamma_rel=1.0+E_0/rest_mass_energy;
    beta_r= v_0/C;
    b_zero=pow(1.0-sqrt(1.0-beta_r*beta_r),2);  // for the Moller cross section
   
    p_0=gamma_rel*Projectile_Mass*v_0;
  //  UMax=p[Param_Offset + 35]; this parameter is not used anymore
    BE_for_exchange=p[Param_Offset + 36];
    BE_for_exchange=BE_for_exchange/Hartree;
    int ExchangeCorrectionMode= p[Param_Offset + 37];
    if(ExchangeCorrection >0)  ExchangeCorrection+=ExchangeCorrectionMode;//ExchangeCorrection=0 no exchange cor.
                                                                          //ExchangeCorrection=1 Ashley method
                                                                          //ExchangeCorrection=2 SBethe method
    DebugMode= p[Param_Offset + 38];
    if(DebugMode )printf("debug mode is on\n");
	return 0;//end copyP_to_Vars
}

    

double DIIMFP_at_omega(double omega)
{
    w_global = omega;
    double q2used, Qrecoil_min, Qrecoil_max;
    
    bool neutral=false;
    double E_1=E_0-omega;
    if(E_1 <=0.0) return 0.0;
 
    //double total_rel_energy= (E_1+rest_mass_energy); //corresponding to E_1
    double rel_cor_factor_1 =  (1 + E_1 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_1 / (Projectile_Mass*C*C)),2);	

    double v_1 = sqrt(2.0*  rel_cor_factor_1*E_1/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double gamma_relativistic_1=1.0+E_1/rest_mass_energy;
    double p_1=gamma_relativistic_1*Projectile_Mass*v_1;
   
    double q1 = p_0 - p_1; //integration boundaries
    double q2 = p_0 + p_1;
    
    Qrecoil_min = recoil_energy_free_electron(q1);
    if (Qrecoil_min > E_0)
    {  
        return 0.0;
    }
    Qrecoil_max= recoil_energy_free_electron(q2);
   

    if (Qrecoil_max > 2*omega) //calculate the momentum of an electron with energy omega, i.e. the maximum transferred moentum
    {  
         double totalE=omega+C*C;
         q2 =sqrt(pow(totalE,2)-pow(C,4))/C+5.0;// the additional amount 5 is because electrons are not stationary, so better go a bit further
       
    }
    if (q2 > q1* q2maxFactor)
        q2used = q1* q2maxFactor;
    else
        q2used = q2; 
    if(q1 <  q_lower) q1=q_lower;  //this is for when we want to calculate partial diimfp's
    if(q2used <   q_lower) q2used=q_lower;
    if(q1 >  q_upper) q1=q_upper;
    if(q2used > q_upper)  q2used=q_upper;
    if (errno !=0) my_perror("DIIMFP_at_omega: before quanc8 an error occured");
    if (neutral && (Projectile_Mass > 1.0) )
        quanc8(fun_neutral, q1, q2used, abserr, relerr);
    //else if(ExchangeCorrection>0 )
        //quanc8(mylossfun_exchange, q1, q2used, abserr, relerr);
    else
        quanc8(mylossfun, q1, q2used, abserr, relerr); 
    double DIIMFP_au=2.0* quanc8result /(pi*v_0*v_0);
    if(DebugMode)
    {
        if(quanc8flag > 0.0)
        {
            printf("omega %6.3f, quanc8flag %9.6f q1 %6.4g q2 %6.4g\n",omega*Hartree,quanc8flag,q1,q2used);
        } 
        if (errno !=0) printf("omega %6.4f\n",omega);
        if (errno !=0) my_perror("DIIMFP_at_omega: an error occured");
    }   
    if(DIIMFP_au < 1e-200) DIIMFP_au=1e-200;
    
    
    return DIIMFP_au;
}

 
double  calc_DSEP(double *DSEPresult, double BeamE,  double theta )
{  
	double  omega, q2used,sep,current_stepsize; 
	double prefactor;
	theta_global=theta;

	omega=FirstEnergy;
	double LastEnergy=FirstEnergy+NStep*StepSize;
	prefactor=2/(pi*v_0*v_0*cos(theta));
	current_stepsize=StepSize;
	if(LastEnergy > 0.9*BeamE)LastEnergy=0.9*BeamE;
	sep=0.0;
	for(int i = 0; omega<LastEnergy; i++)
	{   
		double E_1=E_0-omega;
		w_global = omega;
		double rel_cor_factor_1 =  (1 + E_1 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_1 / (Projectile_Mass*C*C)),2);	

		double v_1 = sqrt(2.0*  rel_cor_factor_1*E_1/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
		double gamma_relativistic_1=1.0+E_1/rest_mass_energy;
		double p_1=gamma_relativistic_1*Projectile_Mass*v_1;
		double q1 = p_0 -p_1; 

		double q2 = p_0 + p_1;
	  
		if (q2 > q1* q2maxFactor)
			q2used = q1* q2maxFactor;
		else
			q2used = q2;
		if(Local_Desep)
		{   
			DSEPresult[i]=0.0;
			for (I_Os_global = 0; I_Os_global< MAXOSC;I_Os_global++)
			{  
				if(abs(Ai[I_Os_global]) > 1e-50)
				{
				   quanc8(DSEPFun_per_Os, q1, q2used, abserr, relerr);
				  
					DSEPresult[i] += Ai[I_Os_global]*prefactor*quanc8result/ Hartree;  // now in eV-1
				}
			}
		}
		else
		{   
			quanc8(DSEPfun, q1, q2used, abserr, relerr);
			
			DSEPresult[i] = prefactor* quanc8result/ Hartree;  // now in eV-1
			
		}
		
		sep=sep+DSEPresult[i]*current_stepsize*Hartree;
		omega=omega+current_stepsize;
				
	}//end loop over energy losses
	return sep;
}
double DDCS_excl_retardation(double omega, double theta)
{   
 
    dcomp eps, oneovereps;  
    double  E_1= E_0-omega;
    double rel_cor_factor_1 =  (1.0 + E_1 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_1 / (Projectile_Mass*C*C)),2);	

    double v_1 = sqrt(2.0*  rel_cor_factor_1*E_1/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double gamma_relativistic_1=1.0+E_1/rest_mass_energy;
    double p_1=gamma_relativistic_1*Projectile_Mass*v_1;
    double q=sqrt(p_1*p_1+p_0*p_0- 2*p_0*p_1*cos(theta));
     if(q/p_0 < 1e-6)  // small angle limit, maybe previous q not accurate due to small difference large numbers (seems not to be an issue)
    {
        
        double q_along = (p_0-p_1*cos(theta));
        double q_perp= p_1*sin(theta);
        double q_small_angle=sqrt(q_perp*q_perp+q_along*q_along);
        q=q_small_angle;
    }
    eps= calc_total_eps( omega, q);
    oneovereps=dcomp(1.0,0.0)/ eps;
    double ddcs=-oneovereps.imag()*p_0*p_0/(q*q);
    return ddcs;  // note  prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity) is added at function call;
}

double DDCS_incl_retardation(double omega, double theta)
{   
    double  v_over_c_square= v_0*v_0/(C*C);
    dcomp eps, oneovereps;  
    double  E_1= E_0-omega;
    double rel_cor_factor_1 =  (1.0 + E_1 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_1 / (Projectile_Mass*C*C)),2);	

    double v_1 = sqrt(2.0*  rel_cor_factor_1*E_1/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double gamma_relativistic_1=1.0+E_1/rest_mass_energy;
    double p_1=gamma_relativistic_1*Projectile_Mass*v_1;
    double theta_e= (p_0-p_1)/p_0; 

    double q=sqrt(p_1*p_1+p_0*p_0- 2*p_0*p_1*cos(theta));
    if(q/p_0 < 1e-6)  // small angle limit, maybe previous q not accurate due to small difference large numbers (seems not to be an issue)
    {
        
        double q_along = (p_0-p_1*cos(theta));
        double q_perp= p_1*sin(theta);
        double q_small_angle=sqrt(q_perp*q_perp+q_along*q_along);
        q=q_small_angle;
    }
    eps= calc_total_eps( omega, q);
    oneovereps=dcomp(1.0,0.0)/ eps;
    double a = v_over_c_square*eps.real()-1.0;
    double b=eps.imag()*v_over_c_square;
    double d=theta*theta-theta_e*theta_e * a ;
    double e= theta_e*theta_e*eps.imag()*v_over_c_square;
    double top= theta*theta+theta_e*theta_e*(a*a+b*b);
    double bottom=d*d+e*e;
    double first= -oneovereps.imag();  //Egerton 3.70
    double ddcs=first*top/bottom;
    return ddcs;  // note  prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity) is added at function call;
}


     

extern "C"
{
	int  Eps1Eps2(double *p, double *Eps1, double *Eps2, double q, int modelchoice)
    {
        dcomp totaleps;
		copyP_to_Vars(p, modelchoice);

		for (int i = 0; i < NStep; i++)
		{  
			double omega = FirstEnergy + StepSize*i;
            totaleps=calc_total_eps( omega, q);

            Eps1[i] = totaleps.real();
            Eps2[i] = totaleps.imag();
		}
       
        my_perror("eps1eps2 an error occured");
		return 0;
	}
    
    int  Eps1Eps2_q(double *p, double *Eps1, double *Eps2, double omega, int modelchoice)

	{
		dcomp eps;
       
        omega=omega/Hartree;
       
		copyP_to_Vars(p, modelchoice);  
        NStep=LastMomentum/Stepsize_qplot;
		for (int i = 0; i < NStep; i++)
		{
			double q =  Stepsize_qplot*i+0.01;  //q= 0 has sometimes troubles
			eps=calc_total_eps( omega, q);//calculated at omega_global
			Eps1[i] = eps.real();
			Eps2[i] = eps.imag();
		}
        if (errno !=0) my_perror(" Eps1Eps2_q an error occured");
		return 0;
	}
    
    int Loss_wide(double *p, double *w_array, double *ELF,  double q, int modelchoice, double *pResultArray)
    {
		//also calculates the constants for the simple DL formulae
        dcomp eps;
        double sumr1 = 0.0;
        double sumr2 = 0.0;
        double sumr3 = 0.0;
        double sumr4 = 0.0;
        double sumr5 = 0.0;
        double sumr6 = 0.0;
        double sumBethe=0.0;
       
		copyP_to_Vars(p, modelchoice);  
        double omega = -0.125/Hartree;  // working in atomic units
        double step = 0.25/Hartree;
       
		for (int i = 0; i < 1300; i++)
		{   
			omega += step;
            eps=calc_total_eps( omega, q);  
            ELF[i]= eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
            w_array[i] = omega*Hartree;  // output in eV
            sumr1 += 2.0 / pi * ELF[i] *  log(omega) * step; 
            sumr2 += 2.0 / pi * ELF[i] * step;
            sumr3 += 2.0 / pi * ELF[i] * omega * log(omega) * step; 
            sumr4 += 2.0 / pi * ELF[i] * omega * step;
            sumr5 += 2.0 / pi * ELF[i] * omega * omega * log(omega) * step; 
            sumr6 += 2.0 / pi * ELF[i] * omega * omega * step;
            sumBethe += 1.0/ (2.0*pi*pi)*ELF[i] * omega * step;
            
            step = step * 1.005;
		}
        pResultArray[0] = sumr2 / exp(sumr1 / sumr2);  // self.C0
        pResultArray[1] = sumr4 / pow(exp(sumr3 / sumr4),2);//self.C1
        pResultArray[2] = sumr6 / pow(exp(sumr5 / sumr6),3); // self.C2
        pResultArray[3] = Hartree * exp(sumr1 / sumr2);  //self.I0 in eV
        pResultArray[4] = Hartree * exp(sumr3 / sumr4);  //self.MIE in eV
        pResultArray[5] = Hartree * exp(sumr4 / sumr5);  //self.Istraggling in eV
        pResultArray[6] =  sumBethe;
        if (errno !=0) my_perror("Loss_wide an error occured");
        return 0;
	}
	
    
    
    int  SurfLossFunc(double *p, double *SurfLoss, double q, int modelchoice)
//gos contributions are not included
	{   
	
		dcomp eps, epsfraction;

		copyP_to_Vars(p, modelchoice);  
		if(Local_Desep)// in this case add the surf loss function of each component weighted by amplitude
		//problems for Drude
		{ 
			for (int i = 0; i < NStep; i++)
			{  
				SurfLoss[i]=0.0;
				w_global = FirstEnergy + StepSize*i;
				
				for(int j=0; j <MAXOSC ; j++)
				{   
					if(abs(Ai[j])>1e-50)
					{
						eps =  get_eps_component(j,q);
						epsfraction = (eps - dcomp(1, 0))*(eps - dcomp(1, 0)) / (eps*(eps + dcomp(1, 0)));
						
						SurfLoss[i]+=Ai[j]*epsfraction.imag();
						
					}
				}
			}
		}
		else // here we calculate the eps of all oscillators combined and get the surf. excitation from that
		{ 
			for (int i = 0; i < NStep; i++)
			{
				double omega = FirstEnergy + StepSize*i;
				SurfLoss[i]=Surface_loss_only(omega, q);
			}
		}
        if (errno !=0) my_perror("SurfLossFunc: an error occured");
		return 0;
	}
    

    double lossfunction_inclGOS_atE(double *p,  double omega, int modelchoice)
	{
		double  GOS,KanekoLoss, Loss_at_E;
		dcomp eps,oneovereps;
	   
		copyP_to_Vars(p, modelchoice);  
 
        w_global = omega/Hartree;
        eps=calculate_eps_osc(0.01);
        oneovereps=dcomp(1.0,0)/eps;
        Loss_at_E = -oneovereps.imag();
    
        GOS=calculate_loss_GOS_dens(0.01);
   
        Loss_at_E +=GOS;
        KanekoLoss=calculate_Loss_AA_LL(0.01);
        Loss_at_E +=KanekoLoss;
    //    if (errno !=0) my_perror(" lossfunction_inclGOS_atE an error occured");
		return(Loss_at_E);
	}
    
    int  lossfunction_inclGOS(double *p, double *lossfunction, double q, int modelchoice)
	{
		
		dcomp eps,oneovereps, chi_Belk;
		copyP_to_Vars(p, modelchoice);  
 
        for (int i = 0; i < NStep; i++)
		{   
			double omega = FirstEnergy + StepSize*i;
			w_global = omega;
            eps=calculate_eps_osc(q);
            oneovereps=dcomp(1.0,0.0)/eps;
			lossfunction[i] =-oneovereps.imag();

            lossfunction[i] +=calculate_loss_GOS_dens(q);

            chi_Belk=calculate_chi_Belkacem( q);
            oneovereps=dcomp(1.0,0.0)/(chi_Belk+dcomp(1.0,0.0));
            lossfunction[i] -=oneovereps.imag();
            lossfunction[i] +=calculate_Loss_AA_LL(q);
            
        }
       // if (errno !=0) my_perror(" lossfunction_inclGOS an error occured");
		return 0;
	}



   // 
   int  DIIMFP(double *p, double *DIIMFPresult,   int modelchoice, double *pResultArray)

	{
        //diimfp calculated for the energy range as defined in chapidif.
        //if this range is smaller than the range that the DIIMFP is non-zero, 
        //then the obtained IMFP will be too large, and straggling and stopping too small 
        double inv_lambda=1E-99;  // so we never get divide by 0
        double stopping=0.0;
        double straggling=0.0;
        double DIIMFP_au=0.0;
       
 		copyP_to_Vars(p,  modelchoice); 
        for (int i = 0; i < NStep; i++)
		{
			double omega    = FirstEnergy + StepSize*i;
            if(ExchangeCorrection==0)
            {
                 DIIMFP_au= DIIMFP_at_omega(omega);
            }
            else if  (omega > 0.5*(E_0+BE_for_exchange))
            {
                 DIIMFP_au=0.0;
            }
            else if (ExchangeCorrection==1)
            {
                DIIMFP_au= DIIMFP_at_omega(omega);
                double DIIMFP_au_exchange= DIIMFP_at_omega(E_0-omega+BE_for_exchange);
                DIIMFP_au=DIIMFP_au+DIIMFP_au_exchange-sqrt(DIIMFP_au*DIIMFP_au_exchange);
            }
            else if (ExchangeCorrection==2)
            {
                double MoellerFactor= 1.0+ pow((omega/(E_0-omega+BE_for_exchange)),2);
                MoellerFactor += -(1.0-b_zero)*omega/(E_0-omega+BE_for_exchange)+b_zero*omega*omega/(E_0*E_0);
                DIIMFP_au= DIIMFP_at_omega(omega)*MoellerFactor;
            }
            DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);// now eV/angstrom
            inv_lambda     += DIIMFP_au * StepSize;
            stopping       += DIIMFP_au * StepSize * omega;
            straggling     += DIIMFP_au * StepSize * omega * omega;
		}
        // the next quantities are only valid estimates if the energy loss range is large enough 
        pResultArray[0]=(1.0/inv_lambda)*BohrRadius;  //IMFP now in ;
        pResultArray[1]=Hartree* (stopping)/BohrRadius;  //stopping now in eV/Angstrom
        pResultArray[2]= Hartree*Hartree* (straggling)/BohrRadius;  //straggling  now in eV^2/Angstrom
        if (errno !=0) my_perror("DIIMFP: an error occured");
        return 0;
	}
	
    int  DIIMFP_for_stopping(double *p, int modelchoice, double *pResultArray)
	{
	   //diimfp calculated for all energy losses with variable stepsize, for stopping straggling and IMFP calculation.
	  
	   double inv_lambda=1E-99;  // so we never get divide by 0
	   double stopping=0.0;
	   double straggling=0.0;
	   copyP_to_Vars(p,  modelchoice); 
	   double omega=0.5*StepSize;
	   double Current_StepSize=StepSize;
	   double DIIMFP_au;
       double omega_max= 2*beta_r*beta_r*gamma_rel*gamma_rel*C*C;

       double MottFactor= (1- beta_r*beta_r*omega/omega_max); // correction factor SalvatPRA eq.16, not sure about this one, esp. for electrons
       
	   bool KeepGoing=true;
       
	   while (KeepGoing)
	   {     
           
            if (ExchangeCorrection==1)
            {
                DIIMFP_au= DIIMFP_at_omega(omega);
                double DIIMFP_au_exchange= DIIMFP_at_omega(E_0-omega+BE_for_exchange);
                DIIMFP_au=DIIMFP_au+DIIMFP_au_exchange-sqrt(DIIMFP_au*DIIMFP_au_exchange);
            }
            else if (ExchangeCorrection==2)
            {
                    double MoellerFactor= 1.0+ pow((omega/(E_0-omega+BE_for_exchange)),2);
                    MoellerFactor += -(1.0-b_zero)*omega/(E_0-omega+BE_for_exchange)+b_zero*omega*omega/(E_0*E_0);
                    DIIMFP_au= DIIMFP_at_omega(omega)*MoellerFactor;
            }
            else  // no exchange
            {
                 DIIMFP_au= DIIMFP_at_omega(omega);
            }
			if(MottCorrection)DIIMFP_au *= MottFactor; // correction factor Salvat PRA  106 032809 eq.16
			inv_lambda     += DIIMFP_au * Current_StepSize;
			stopping       += DIIMFP_au * Current_StepSize * omega;
			straggling     += DIIMFP_au * Current_StepSize * omega * omega;
            
			Current_StepSize=StepSize+lin_cont_deltaE*omega;
			omega=omega+Current_StepSize;
            
            if (omega < omega_max) MottFactor= (1- beta_r*beta_r*omega/omega_max);// keep Mott factor the same for omega > omega_max
         

			if((ExchangeCorrection  != 0) && (omega > 0.5*(E_0+BE_for_exchange))) KeepGoing=false; 
            if(omega > 1.5* omega_max+2.0) KeepGoing=false; 
		}
        if (errno !=0) my_perror("halfway _for_stopping: an error occured");
		pResultArray[0]=(1.0/(inv_lambda+1.0e-99))*BohrRadius;  //IMFP now in Angstrom;
		pResultArray[1]=Hartree* (stopping)/BohrRadius;  //stopping now in eV/Angstrom
		pResultArray[2]= Hartree*Hartree* (straggling)/BohrRadius;  //straggling  now in eV^2/Angstrom
        if(DebugMode)
        {
            if (errno !=0) my_perror("DIIMFP_for_stopping at end: an error occured");
        }
		return 0;
	}
	
	
	double  DSEP(double *p, double *DSEPresult,  int modelchoice)
	{
        //projectile=1: electron  projectile = 2: proton
		double sep;
 		copyP_to_Vars(p,  modelchoice); 
		sep=calc_DSEP(DSEPresult,E_0,theta0);
        if (errno !=0) my_perror("DSEP: an error occured");
        return sep;
	}
    
	void GOS_Scaling_init(double *p)  
	{   int  modelchoice=1;  // modelchoice does not matter here, but has to be defined
		copyP_to_Vars(p,  modelchoice); 
		double CurrentFSum, CurrentE,Currentq,stepsize_au;
		for (int i = 0; i < MAXGOS; i++)
		{
			if (AiGOS[i] > 0)
			{
				for(int iq=0;iq<80;iq++)
				{   
					Currentq = GOSScalingStepSize*iq;                
					stepsize_au=Edgei_GOS[i]/5000;  // energy   step equal to 1/500 of the edge energy
					CurrentE=Edgei_GOS[i]+0.5*stepsize_au;
					CurrentFSum = 0.0;
					while(CurrentE < (200*Edgei_GOS[i]+Currentq*Currentq/2.0))  // so we integrate up to 200 times the  edge energy+q^2/2
					{
						CurrentFSum += GOSx(n_i_GOS[i], l_i_GOS[i], Zi_GOS[i], Edgei_GOS[i], Currentq, CurrentE)*stepsize_au;
						stepsize_au=stepsize_au*1.001;
						CurrentE=CurrentE+stepsize_au;
					}
				   GOS_ScalingFactor[i][iq] = 1.0 / CurrentFSum;
			      printf(" i  %i,Currentq %6.3f,scalingfactor, GOS_ScalingFactor %6.2f\n",i, Currentq, GOS_ScalingFactor[i][iq] );
			   }
			}		
		}
        if (errno !=0) my_perror("an error occured");
        return;
	}
    double TaucSumRule(double  C, double E0,  double Egap)  // this does NOT call copyP_to_vars
    {   
        double sum= TaucNormalisation(  C, E0,  Egap);
        return sum;
    }
    

     int  Kramers_Kronig_eps1_from_eps2(double FirstE, double DeltaE,double  eps1_infinity ,int length, double const *eps2, double *eps1)  // does not depend on  copyP_to_Vars
     {  
        double omega_prime;
        for (int i = 0; i < length; i++)
        {
           double KKsum = 0.0;
           double omega = i * DeltaE + FirstE;  
           for (int k = 0; k < length; k++)
           {
                omega_prime = k * DeltaE + FirstE;
                if (i != k)
                     KKsum += omega_prime*eps2[k] / (omega_prime * omega_prime - omega * omega) * DeltaE;// Wooten eq. 6.52 
           }                                                                                                               
           eps1[i]  =2.0* KKsum / pi + eps1_infinity; //the factor 2 because we go from 0 to large rather than -large to large
        }
        if (errno !=0) my_perror("Kramers_Kronig_eps1_from_eps2: an error occured");
        return 0;
     }  
     
    int  Kramers_Kronig_eps2_from_eps1(double FirstE, double DeltaE,double  eps1_infinity , int length, double const *eps1, double *eps2) // does not depend on  copyP_to_Vars
     {  
        double omega_prime;
        for (int i = 0; i < length; i++)
        {
           double KKsum = 0.0;
           double omega= i * DeltaE + FirstE;  
           for (int k = 0; k < length; k++)
           {
                omega_prime = k * DeltaE + FirstE;
                if (i != k)
                     KKsum += omega * (eps1[k]- eps1_infinity ) / (omega_prime*omega_prime - omega*omega) * DeltaE;// Wooten eq. 6.52 
           }                                                                                                               
           eps2[i]  = - 2.0* KKsum / pi; //the factor 2 because we go from 0 to large rather than -large to large
        }
        if (errno !=0) my_perror("Kramers_Kronig_eps2_from_eps1: an error occured");
        return 0;
     }  
     
      
double   DDCS_longitudinal(double *p, double *DDCS,  int modelchoice, int integrate_phi)
	{
		double theta_step, omega, stopping, DDCS_au;
        int Nthetastep,ienergy,itheta;
        copyP_to_Vars(p,  modelchoice); 
        Nthetastep=LastMomentum/Stepsize_qplot;// so number of steps here similar number of steps omega-q plot
        theta_step=theta_max/Nthetastep;
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        stopping=0.0;
        for(itheta=0; itheta < Nthetastep; itheta++)
        {
            double theta=(0.5+itheta)*theta_step;
           
            for(ienergy =0;ienergy < NStep; ienergy++)
            {   
                omega= FirstEnergy +ienergy*StepSize;
                DDCS_au= prefactor * DDCS_excl_retardation(omega,theta);
                stopping += 2 * pi * sin(theta) * theta_step * DDCS_au * omega * StepSize; 
                DDCS[itheta*NStep+ienergy] = DDCS_au * (BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;
            }
        }

        stopping *= UnitCellDensity;// now in a.u. energy/ a.u. distance
        stopping *= Hartree / BohrRadius ;  // now in eV/Angstrom
        if (errno !=0) my_perror("DDCS_longitudinal: an error occured"); 
        return stopping;
    }
    
double   DDCS_total(double *p, double *DDCS,   int modelchoice, int integrate_phi)  // including Cerenkov, i.e. retardation effects
	{
		double theta_step, omega, stopping, DDCS_au;
        int Nthetastep,ienergy,itheta;
        copyP_to_Vars(p,  modelchoice); 
    
        stopping=0.0;
        Nthetastep=LastMomentum/Stepsize_qplot;// so number of steps here similar number of steps omega-q plot
        theta_step=theta_max/Nthetastep;
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        for(itheta=0; itheta < Nthetastep; itheta++)
        {
            double theta=(0.5+itheta)*theta_step;
            for(ienergy =0;ienergy < NStep; ienergy++)
            {
                omega= FirstEnergy +ienergy*StepSize;
                DDCS_au= prefactor * DDCS_incl_retardation(omega,theta);
                stopping += 2 * pi * sin(theta) * theta_step * DDCS_au * omega * StepSize; 
                DDCS[itheta*NStep+ienergy] = DDCS_au * (BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;
            }
        }
        stopping *=  UnitCellDensity;// now in a.u. energy/ a.u. distance
        stopping *= Hartree / BohrRadius ;  // now in eV/Angstrom
        if (errno !=0) my_perror("DDCS_total: an error occured");    
        return stopping;
    }    
    
int  DDCS_at_theta(double *p, double *DDCS, double *DDCS_incl_ret,  int modelchoice, double theta)  
    {   double  result_au, result_au_incl_ret;
        copyP_to_Vars(p,  modelchoice);
        theta=theta/1000.0; // now in rad  
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
       
      
        for(int ienergy =0;ienergy < NStep; ienergy++)
            {                   
                double omega = FirstEnergy +ienergy*StepSize;
                result_au = prefactor *  DDCS_excl_retardation(omega, theta);  //Egerton 3.32
                DDCS[ienergy]=result_au*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
                result_au_incl_ret = prefactor *  DDCS_incl_retardation(omega, theta);  //Egerton 3.32
                DDCS_incl_ret[ienergy]=result_au_incl_ret*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
            }    
        if (errno !=0) my_perror("DDCS_at_theta: an error occured");
        return 0;    
    } 
    
int  DDCS_at_omega(double *p, double *xaxis, double *DDCS,  double *DDCS_incl_ret,  int modelchoice,  double omega)
	{
		double  result_au, result_au_incl_ret;
        int Nthetastep,itheta;
        copyP_to_Vars(p,  modelchoice); 
        Nthetastep=LastMomentum/Stepsize_qplot;// so number of steps here is the same as the number of steps omega-q plot
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        omega=omega/Hartree;
        for(itheta=0; itheta < Nthetastep; itheta++)
        {
            double theta=xaxis[itheta]/1000.0;   // in rad          
            
            result_au = prefactor *  DDCS_excl_retardation(omega, theta);  //Egerton 3.32
            DDCS[itheta]=result_au*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
            result_au_incl_ret = prefactor *  DDCS_incl_retardation(omega, theta);  //Egerton 3.32
            DDCS_incl_ret[itheta]=result_au_incl_ret*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
        }
        if (errno !=0) my_perror("DDCS_at_omega: an error occured");
        return 0;
    }
    int  DCS(double *p, double *xaxis, double *DCS,  double *DCS_incl_ret,    int modelchoice)
	{
	
        
        int Nthetastep,itheta;
        copyP_to_Vars(p,  modelchoice); 
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        
        Nthetastep=LastMomentum/Stepsize_qplot;// so number of steps here is the same as the number of steps omega-q plot
        for(itheta=0; itheta < Nthetastep; itheta++)
        {
             DCS[itheta]=0;
             DCS_incl_ret[itheta]=0;
            
        
            for(int ienergy =0;ienergy < NStep; ienergy++)    
            {    
                double omega= FirstEnergy +ienergy*StepSize;
                double theta=xaxis[itheta]/1000.0;   // in rad  
                double result_au = prefactor *  DDCS_excl_retardation(omega, theta);
                DCS[itheta] += result_au * StepSize_eV*(BohrRadius*BohrRadius)/Hartree;
                result_au = prefactor *  DDCS_incl_retardation(omega, theta);
                DCS_incl_ret[itheta] += result_au * StepSize_eV*(BohrRadius*BohrRadius)/Hartree;
            }
        }
        if (errno !=0) my_perror("DCS: an error occured");
        return 0;
    }
    
    
     int  calc_REELS(double *p,   double StartReels, int N_REELS_Step, int modelchoice,  double *REELSresult, int EELS, double EELSPathLength)
//works in eV
	{
        
	double lambda; 
	//projectile=1: electron  projectile = 2: proton
	double GE,arg;
    double sprob[3];

	copyP_to_Vars(p,  modelchoice); 
    double  *DIIMFPresult = new double [NStep];
    double  *DSEPresultIn = new double [NStep];
    double  *DSEPresultOut = new double [NStep];
    double  *NormDSEPresult = new double [NStep];
    double  *normDIIMFP = new double[NStep];
    double  *vSignalconv = new double[2*NStep];

      // Dimensions of the 3D array    
    int MaxPartInt=10;
    double *PartialIntensity= new double [MaxPartInt];
    int MaxSurfEx=3;

   
    double*** matrix2 = new double**[N_REELS_Step];
 
    for (int i = 0; i < N_REELS_Step; i++) {
 
        // Allocate memory blocks for
        // rows of each 2D array
        matrix2[i] = new double*[MaxPartInt];
 
        for (int j = 0; j < MaxPartInt; j++) {
 
            // Allocate memory blocks for
            // columns of each 2D array
            matrix2[i][j] = new double[MaxSurfEx];
        }
    }
    FirstEnergy=0.5* StepSize;  // This makes sure FirstEnergy is sensible, not what in user interface
    double inv_lambda  =1.0e-99;
    for (int i = 0; i < NStep; i++)
		{
			double omega    = FirstEnergy + StepSize*i;
            double DIIMFP_au= DIIMFP_at_omega(omega);
            DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);// now eV/angstrom
            inv_lambda     +=  DIIMFPresult[i]*StepSize_eV;// now in angstrom^{-1}
        }
    lambda=1.0/inv_lambda;
    if(EELS == 0)
    {   
        for (int j = 0; j < MaxPartInt; j++)
            PartialIntensity[j] = 1.0 + PIcoef1*j + PIcoef2*j*j + PIcoef3*j*j*j; 
    }   
    else
    {   
        double Lambda_Poisson=EELSPathLength/lambda;
        double zero_loss=Poisson(0, Lambda_Poisson);
        for (int j = 0; j < MaxPartInt; j++)
        {   
            PartialIntensity[j] =Poisson(j, Lambda_Poisson)/zero_loss;  // so elastic peak intensity = 1 
        }   
    } 
    for (int i=0; i< NStep; i++)
    {
        normDIIMFP[i]=DIIMFPresult[i]*lambda*fraction_DIIMFP;
    }
    double SurfExProbin  = calc_DSEP(DSEPresultIn,E_0, theta0);
    double SurfExProbout = calc_DSEP(DSEPresultOut,E_0, theta1);
    if(SurfExProbin > 1.0) SurfExProbin=1.0;
    if(SurfExProbout> 1.0) SurfExProbout=1.0;

    if(surf_ex_factor > 0.0 && SurfExProbin!=0.0 && SurfExProbout !=0.0)
    {
        for (int i=0; i< NStep; i++)
        { // take the average of the two similar DSEPs
            NormDSEPresult[i]=0.5*StepSize*Hartree*(DSEPresultIn[i]/SurfExProbin+DSEPresultOut[i]/SurfExProbout);
        }
        SurfExProbin=SurfExProbin*surf_ex_factor;
        SurfExProbout=SurfExProbout*surf_ex_factor;
        sprob[0]= (1.0-SurfExProbin)*(1.0-SurfExProbout);
        sprob[1]= (1.0-SurfExProbin)*SurfExProbout+(SurfExProbin)*(1.0-SurfExProbout);
        sprob[2]= SurfExProbin*SurfExProbout;
    }
     else
     {
        for (int i=0; i< NStep; i++)
        { // take the average of the two similar DSEPs
            NormDSEPresult[i]=0.0;
        }
        sprob[0]=1.0; sprob[1]=0.0; sprob[2]=0.0;
    }
   

    for (int ipartint = 0; ipartint <MaxPartInt; ipartint++)
	{
		if (ipartint == 0)  // fill up the first element with the spectrometer response,  a Gaussian
		{
			for (int i = 0; i < N_REELS_Step; i++)
			{
				GE = StartReels + StepSize_eV*i;

				//GE = StepSize_eV*(i - 20);
				arg = GE*GE / (2.0*sigma*sigma);
				if (arg < 40)
				{
					matrix2[i][0][0] = exp(-arg); // take maximum to be 1, norm exp elastic peak to 1 for easy comp
				}
				else
				{
					matrix2[i][0][0] = 0.0;
				}
			}
			
		}
		else  //  fill up subsequent elements with  a Gaussian convoluted i times with NDIIMFP
		{
			for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
			for (int i = 0; i < NStep; i++)
				for (int j = 0; j < NStep; j++){
					vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint - 1][0] * normDIIMFP[j] * StepSize_eV;
				}

			// put in matrix2[i][ipartint][0]
			for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][0] = vSignalconv[j];
		}
// add convolution with normalised DSEP to next columns    

        for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
        for (int i = 0; i < NStep; i++){
            
            for (int j = 0; j < NStep; j++){
                vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint][0] * NormDSEPresult[j];
            }
        }
        // put in matrix2[i][ipartint][1]
        for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][1] = vSignalconv[j];


        for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
        for (int i = 0; i < NStep; i++){
            for (int j = 0; j < NStep; j++){
                vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint][1] * NormDSEPresult[j];
            }
        }
        // put in matrix2[i][ipartint][2]
        for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][2] = vSignalconv[j];
    }
    fstream myfile;
    myfile.open("loss_dist.txt",fstream::out);
    myfile << "Energy:\t No surf loss, \t surf loss\t two surf loss" <<std::endl; 
    
    for (int j = 0; j < N_REELS_Step; j++)
	{
        double no_surfplasmon,onesurfplasmon,twosurfplasmon; 
        no_surfplasmon=0.0;onesurfplasmon=0.0;twosurfplasmon=0.0; 
		for (int jj = 0; jj < MaxPartInt; jj++)
        {
            no_surfplasmon+=PartialIntensity[jj] * matrix2[j][jj][0]; // for convenience we keep max height elastic peak at 1
            onesurfplasmon+=(sprob[1]/sprob[0])*PartialIntensity[jj] * matrix2[j][jj][1];
            twosurfplasmon+= (sprob[2]/sprob[0])*PartialIntensity[jj] * matrix2[j][jj][2];
        }
        REELSresult[j] =  no_surfplasmon+onesurfplasmon+twosurfplasmon; 
        myfile << StartReels+j*StepSize_eV << "\t"<< no_surfplasmon << "\t"<<onesurfplasmon<< "\t"<<twosurfplasmon<<std::endl;  
    }    
    myfile.close();
    myfile.open("PartInt.txt",fstream::out);
    myfile << "Part.Int.:\t"; 
    for (int j=0; j< MaxPartInt;j++)// Prints row of x
    {        
        myfile << j << "\t";  
    }

    myfile<< std::endl;

    for (int j=0; j< N_REELS_Step;j++) //This variable is for each row below the x 
    {        
        myfile << StartReels+j*StepSize_eV << "\t";

        for (int jj=0; jj<MaxPartInt;jj++)
        {                      
            myfile <<matrix2[j][jj][0] << "\t";
        }
        myfile<<std::endl;
    }
    myfile.close();

    delete[] NormDSEPresult;
    delete[] DSEPresultIn;
    delete[] DSEPresultOut;
    delete[] DIIMFPresult;
    delete[] vSignalconv;
    delete[] normDIIMFP;
    for (int i = 0; i < N_REELS_Step; i++) {
        for (int j = 0; j < MaxPartInt; j++) {
            delete[] matrix2[i][j];
        }
        delete[] matrix2[i];
    }
    delete[] matrix2;
    if (errno !=0) my_perror("calc_REELS: an error occured");
    return 0;
    }
	 
}  //matching brace of extern C
