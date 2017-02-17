#load and format data and parameter sets from Teixeira paper
library(R.matlab)
library(deSolve)

#load and reshape parameter estimates etc from Teixeira model
ode_pars <- readMat('~/../Dropbox/USFpostdoc/WALB_DEB/data/ode_pars.mat')[[1]]
names(ode_pars) <- sub("\\.", "_", dimnames(ode_pars)[[1]])
times <- as.vector(ode_pars['a'][[1]])
pars <- unlist(ode_pars[c(6:16)])
inits <- unlist(ode_pars[c(2:5)])
states <- c(H = unname(inits['E_Hh']), E = unname(inits['E_h']), L=unname(inits['L_h']), f_n=(-1e-10 * unname(pars['f_slope']) + unname(pars['f_intercept'])))
states_log <- c(H = unname(inits['E_Hh']), E = unname(inits['E_h']), L=unname(inits['L_h']), f_n=1.696393085887949)
params = pars
params['f_slope'] <- -params['f_slope']

#load data
obs_table <- read.csv('~/bitbucket/walb_deb/data/teixeira2014tab02.csv', skip=4)
#rename to match model output
names(obs_table) <- c("time", "Ww", "L_cul", "TBW", "wdratio")
#change units to g
obs_table$Ww <- obs_table$Ww*1000
#reshape and reorder columns
obs_long <- reshape2::melt(obs_table, id="time")[,c(2,1,3)]
