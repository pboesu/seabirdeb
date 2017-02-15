/* file teix_simple.c */
#include <R.h>
#include <math.h>

static double parms[10];
#define L_m   parms[0]
#define p_Am  parms[1]
#define v     parms[2]
#define k_J   parms[3]
#define kap   parms[4]
#define T_A   parms[5]
#define T_ref parms[6]
#define T_b   parms[7]
#define E_G   parms[8]
#define f_slope parms[9]
/* #define f_intercept parms[10] */


/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=10;
  odeparms(&N, parms);
}

/* Derivatives and 2 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  if (ip[0] <1) error("nout should be at least 1");
    double E_m;
    double g;
    double TC;
    double vT;
    double kT_J;
    double pT_Am;
    double e;
    double rT;
    double pT_C;
    #define H   y[0]
    #define E   y[1]
    #define L   y[2]
    #define f_n y[3]
    #define wdratio y[4]

    E_m = p_Am/ v; //#% J/cm^3, reserve capacity [E_m]
    g = E_G/ kap/ E_m; //#% -, energy investment ratio
    TC = exp(T_A/T_ref - T_A/T_b);
    vT = v * TC;
    kT_J = k_J * TC;
    pT_Am = p_Am * TC;
    e = vT * E/ pow(L, 3)/ pT_Am; //#% -, scaled reserve density;
    rT = vT * (e/ L - 1/ L_m)/ (e + g); //#% 1/d, spec growth rate
    pT_C = E * (vT/ L - rT); //#% J/d, scaled mobilisation

    /* derivatives */
    ydot[0] = (1 - kap) * pT_C - kT_J * H; // dH J
    ydot[1] = pT_Am * f_n * pow(L, 2) - pT_C; //dE
    ydot[2] = rT * L/3; //dL
    ydot[3] = -f_slope; //df_n

    /* limit the functional response to >= 0 */

      if(y[3] < 0)
           ydot[3] = 0;

    /* calculate derived quantities: culmen length and dry weight */
      yout[0] = y[2]*1.074;

    //# #wet weight. this is all hard coded now, should not be!
      static const double w_E = 23.9; // # molecular weight of reserve g mol^-1
      static const double d_v = 0.5; // # specific density of structure
      static const double mu_E = 550000; // # chemical potential of reserve J / mol
      /* since I can't quite figure out the way to multiply in time-dependent wdratios,
       * we are sticking to dry weight instead
       *
       * IDEA: make wdratio a state variable, then the multiplication should simply work
       * using y[4]
       */
      //double wdratio = -1.37e-3 * *t + 2.09;
      double omega = p_Am * w_E / (v * d_v * mu_E);
      yout[1] = pow(y[2], 3) * (1 + y[4] * omega) * d_v;
}

/* The Jacobian matrix */
/* void jac(int *neq, double *t, double *y, int *ml, int *mu,
          double *pd, int *nrowpd, double *yout, int *ip)
 {
   pd[0]               = -k1;
   pd[1]               = k1;
   pd[2]               = 0.0;
   pd[(*nrowpd)]       = k2*y[2];
   pd[(*nrowpd) + 1]   = -k2*y[2] - 2*k3*y[1];
   pd[(*nrowpd) + 2]   = 2*k3*y[1];
   pd[(*nrowpd)*2]     = k2*y[1];
   pd[2*(*nrowpd) + 1] = -k2 * y[1];
   pd[2*(*nrowpd) + 2] = 0.0;
 } */
/* END file mymod.c */
