#' Teixeira albatross model
#'
#' @param t vector of times
#' @param states vector of initial state values
#' @param params vector of parameters
#'
#' @return derivatives of the ODE
#' @export
#'
walb_deb <- function(t, states, params){
  with(as.list(c(states, params)), {

    E_m = p_Am/ v;     #% J/cm^3, reserve capacity [E_m]
    g = E_G/ kap/ E_m; #% -, energy investment ratio

    TC = exp(T_A/T_ref - T_A/T_b);
    vT = v * TC;
    kT_J = k_J * TC;
    pT_Am = p_Am * TC;

    e = vT * E/ L^3/ pT_Am;             #% -, scaled reserve density;
    rT = vT * (e/ L - 1/ L_m)/ (e + g); #% 1/d, spec growth rate
    pT_C = E * (vT/ L - rT);            #% J/d, scaled mobilisation

    #[f_n] = get_f(t, f_slope, f_intercept);


    #derivatives
    dH = (1 - kap) * pT_C - kT_J * H; #% J
    dE = pT_Am * f_n * L^2 - pT_C;    #% J
    dL = rT * L/3;                    #% cm
    df_n = -f_slope               #% -  #different from matlab implementation

    #limit the functional response to >= 0
    if(states[4] < 0) df_n <- 0

    #return derivatives
    list(c(dH, dE, dL, df_n))
  })
}

#' Teixeira albatross model with logistic food in nestling stage
#'
#' @param t vector of times
#' @param states vector of initial state values
#' @param params vector of parameters
#'
#' @return derivatives of the ODE
#' @export
#'
walb_deb_log_food <- function(t, states, params){
  with(as.list(c(states, params)), {

    E_m = p_Am/ v;     #% J/cm^3, reserve capacity [E_m]
    g = E_G/ kap/ E_m; #% -, energy investment ratio

    TC = exp(T_A/T_ref - T_A/T_b);
    vT = v * TC;
    kT_J = k_J * TC;
    pT_Am = p_Am * TC;

    e = vT * E/ L^3/ pT_Am;             #% -, scaled reserve density;
    rT = vT * (e/ L - 1/ L_m)/ (e + g); #% 1/d, spec growth rate
    pT_C = E * (vT/ L - rT);            #% J/d, scaled mobilisation

    #[f_n] = get_f(t, f_slope, f_intercept);
    log_df_n <- function(t, lAsym=0.334, uAsym=1.75, rate = -0.038, xmid= 85.14){
      rate*exp(-rate*(t - xmid))*(uAsym - lAsym)/(1 + exp( -rate*(t - xmid)))^2
    }


    #derivatives
    dH = (1 - kap) * pT_C - kT_J * H; #% J
    dE = pT_Am * f_n * L^2 - pT_C;    #% J
    dL = rT * L/3;                    #% cm
    df_n = log_df_n(t)               #% -  #different from matlab implementation

    #return derivatives
    list(c(dH, dE, dL, df_n))
  })
}
