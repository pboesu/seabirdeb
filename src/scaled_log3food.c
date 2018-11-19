/* file scaled_std_deb.c */
#include <R.h>
#include <math.h>

static double scaled_std_parms[16];
#define L_m   scaled_std_parms[0]
#define p_Am  scaled_std_parms[1]
#define v     scaled_std_parms[2]
#define k_J   scaled_std_parms[3]
#define kap   scaled_std_parms[4]
#define T_A   scaled_std_parms[5]
#define T_ref scaled_std_parms[6]
#define T_b   scaled_std_parms[7]
#define E_G   scaled_std_parms[8]
#define f_lower scaled_std_parms[9]
#define f_rate scaled_std_parms[10]
#define f_mid scaled_std_parms[11]
#define delta_M scaled_std_parms[12]
#define wdratio_slope scaled_std_parms[13]
#define E_Hb  scaled_std_parms[14]
#define E_Hp  scaled_std_parms[15]
/* #define f_intercept parms[10] */


/* initializer  */
void init_scaled_log3food(void (* odeparms)(int *, double *))
{
  int N=16;
  odeparms(&N, scaled_std_parms);
}

/* Derivatives and 2 output variable */
void d_scaled_log3food(int *neq, double *t, double *y, double *ydot,
                       double *yout, int *ip)
{
  if (ip[0] <1) error("nout should be at least 1");
  double E_m;
  double g;
#define e   y[0]
#define l   y[1]
#define uH   y[2]
#define uR y[3]
#define f y[4]

  E_m = p_Am/ v;
  g = E_G/ kap/ E_m;
  double k_M = v / (g * L_m);
  double l_T = 0;
  double k = k_J / k_M;
  double V_m =pow(v/(k_M*g), 3);

  double u_Hp = E_Hp/(g*E_m*V_m);

  static const double w_E = 23.9; //# molecular weight of reserve g mol^-1
  static const double d_v = 0.16; //# specific density of structure
  static const double mu_E = 550000; //# chemical potential of reserve J / mol
  double omega = p_Am * w_E / (v *d_v * mu_E); //#omega

  const double uAsym = 1;

  /* derivatives */
  ydot[0] = k_M * g * (f - y[0]) / y[0]; //d/dt e
  ydot[1] = k_M / 3 * (y[0] - y[1] - l_T) / (1 + y[0] / g); //d/dt l
  ydot[2] = k_M * (y[2] < u_Hp) * ((1 - kap) * y[0] * pow(y[1], 2) * (g + y[1] + l_T) / (g + y[0]) - k * y[2]); //d/dt uH
  ydot[3] = k_M * (y[2] > u_Hp) * ((1 - kap) * y[0] * pow(y[1], 2) * (g + y[1] + l_T) / (g + y[0]) - k * u_Hp); //d/dt uR
  ydot[4] = f_rate * exp(-f_rate * ( *t - f_mid)) * (uAsym - f_lower)/ pow((1 + exp( -f_rate * ( *t - f_mid))), 2); //d/dt f_n
  ydot[5] = wdratio_slope; //-1.37e-3; //d/dt wdratio;


  /* calculate derived quantities: */
  yout[0] = y[2]*L_m/delta_M; //L_cul culmen length
  yout[1] = pow(y[1]*L_m, 3) * (1 + y[0] * omega) * d_v *y[5]; // Ww wet weight


}
