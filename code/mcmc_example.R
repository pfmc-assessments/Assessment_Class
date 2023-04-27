# Must put ss3 executable in "path" directory with other model files
# There is no way around this given how adnuts is currently written

# In addition to standard input files, you must run the model in MLE mode
# within the "path" directory and have a positive definite hessian

# I *think* the hessian is only necessary for the adnuts wrapper functions,
# since many of the diagnostic functions compare MLE and MCMC output.
# You *should* be able to run MCMC from the command line without a hessian.

# pak::pkg_install('Cole-Monnahan-NOAA/adnuts')

fit <- adnuts::sample_rwm(model = 'ss_win', # this is the name of the executable
                          path =  # directory with executable, input file, MLE output files (including covariance)
                          iter = 200000,
                          thin = 100, # thin to save memory, could try not
                          chains = 3)
# By default the first 50% is burn in.

adnuts::plot_marginals(fit)
summary(fit)
adnuts::launch_shinyadmb(fit)
adnuts::pairs_admb(fit, pars=1:6, order='slow')
adnuts::pairs_admb(fit, pars=1:6, order='slow', diag='hist')
adnuts::plot_uncertainties(fit)


post <- adnuts::extract_samples(fit)
