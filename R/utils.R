################################################################################
########################  Function in senerio summary  #########################
################################################################################


pDLT <- function(proba, wm){

  # Compute the probability of observing a DLT given Probability array
  #
  # \code{pDLT} Compute the probability of observing a DLT Given Probability
  #
  # This is called in \code{\link{nTTP_summary}} for Compute the probability of
  # observing a DLT Given Probability vector
  #
  # @param proba   (numerical vector) The probability vector with the length to
  # be the same and the number of elements in wm
  # @param wm      (numerical matric) Toxicity weighted matrix, with row be
  # the type of the toxicity and column be the toxicity grade
  #
  # @return The probability of observing a DLT given the scenerio
  #
  # @importFrom arrayhelpers vec2array
  # @return the probabililty of observing a DLT

  n.tox.grade <- ncol(wm)  # retrieve the No. of tox grade
  n.tox.type  <- nrow(wm)  # retrieve the NO. of tox type

  DLT_array <- array(NA, c(rep(n.tox.grade, n.tox.type)))
  capacity <- n.tox.grade^n.tox.type
  for(tt in 1 : capacity) {
    index <- vec2array(tt, dim = c(rep(n.tox.grade, n.tox.type)))
    DLT_array[tt] <- as.numeric(max(wm[t(rbind(1 : n.tox.type, index))]) >= 1)
  }
  return(sum(proba * DLT_array))
}

mnTTP <- function(proba, nTTP.all){

  # Compute mean nTTP given probability array
  #
  # \code{mnTTP} Compute mean nTTP given probability array
  #
  # This is called in \code{\link{nTTP_summary}} for Compute Mean nTTP Given
  # Probability Arrayr
  #
  # @param proba   (numerical vector) The probability vector with the length to
  # be the same and the number of elements in wm
  # @param nTTP.all The output of \code{\link{nTTP.array}}
  #
  # @return The mean nTTP given probability array

  return(sum(proba * nTTP.all))
}

################################################################################
########################### Function in simulation  ############################
################################################################################

gen_nTTP_dlt <- function(dose_cycle_PatID, tox.prob.M, wm, nTTP.all,
                         eff.structure, eff.Sigma, eff.sd_trans,
                         patlist){

  ## Generate nTTP and DLT and efficacy
  ## dose_cycle_PatID: dose, cycle, PatID
  ## tox.prob.M:    tox prob matrix, dim: dose, cycle, type, grade
  ## wm:            weighted matrix
  ## nTTP.all:      output generated from function nTTP.array
  ## eff.structure: parameter setting for generating efficacy data (a 6 by 6 matrix)
  ## eff.Sigma:     correlation matrix for efficacy across cycles
  ## eff.sd_trans:  the parameter controls the variance and skewness of the efficacy data
  ## patlist:       the simulation data, for searching the history efficacy data
  #' @importFrom stats pnorm rnorm


  # sample the grade level for each toxiciy type
  n.tox.grade <- dim(tox.prob.M)[4]
  dose <- as.numeric(dose_cycle_PatID[1])
  cycle <- as.numeric(dose_cycle_PatID[2])
  nttp.indices <- apply(tox.prob.M[dose, cycle, , ], 1,
                        function(o){sample(1:n.tox.grade, 1, prob = o)})

  # generate the corresponding dlt and the nTTP
  y.dlt <- as.numeric(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1)
  y.nTTP <- nTTP.all[matrix(nttp.indices, nrow = 1)]

  ## generate efficacy data
  if(cycle == 1){
    # y.effz record the transfer variable
    y.effz <- rnorm(1, mean = eff.structure[dose, 1],  sd = eff.Sigma[1, 1])
    y.eff <- pnorm(y.effz, mean = 0, sd = eff.sd_trans)
  } else{
    ChoLD <- t(chol(eff.Sigma[1:(cycle - 1), 1:(cycle - 1)]))
    CinvC <- forwardsolve(ChoLD, eff.Sigma[cycle, 1:(cycle - 1)])
    y.effz.old <- patlist$effz[which(patlist$PatID == dose_cycle_PatID[3])]
    y.dose.old <- patlist$dose[which(patlist$PatID == dose_cycle_PatID[3])]

    CinvM <- forwardsolve(ChoLD, y.effz.old -
                            eff.structure[cbind(y.dose.old, 1:(cycle - 1))])
    CondMu <- eff.structure[dose, cycle] + sum(CinvC * CinvM)
    CondV <- eff.Sigma[cycle, cycle] - sum(CinvC^2)

    y.effz <- rnorm(1, mean = CondMu, sd = sqrt(CondV))
    y.eff <- pnorm(y.effz, mean = 0, sd = eff.sd_trans)
  }

  return(c(y.nTTP = y.nTTP, y.dlt = y.dlt, y.eff = y.eff, y.effz = y.effz))
}


################################################################################
######################### Function in model fitting  ###########################
################################################################################

phase1stage1 <- function(patlist, ctrl_param, n.iters = 5000, burn.in = 5000,
                         retrieve_param = c("beta_dose", "beta_other"),
                         inits.list.set = list(), n.chains = 1,
                         dose_flag = 0){

  ## model fitting for stage 1 ##
  ## patlist:        the current patients record
  ## ctrl_param:     the prior setting
  ## n.iters:        number of iterations for the MCMC chain after burn in
  ## n.burn:         number of burn in
  ## retrive_param:  the parameters we want a posterior samples for the further study
  ## inits.list.mod: modification of starting point of parameters, currently use default..
  ## n.chains:       number of MCMC chains in the fitting, default = 1
  ## dose_flag:      when dose_flag equals 1, we need to degenerate design matrix
  ## return: return the posterior samples of the MCMC chain
  #' @importFrom rjags jags.model jags.samples
  #' @importFrom stats runif update
  #' @importFrom utils modifyList

  uniq_ID <- unique(patlist$PatID)
  n.subj <- length(uniq_ID)                       ## number of patients
  n <- length(patlist$PatID)                      ## number of records
  if(dose_flag == 1){
    X_y <- cbind(patlist$cycle, patlist$dose)     ## design matrix
  } else {
    X_y <- cbind(1, patlist$cycle, patlist$dose)  ## design matrix
  }

  W_y <- matrix(0, n, n.subj)
  W_y[cbind(1:n, as.numeric(sapply(patlist$PatID, function(a) {
    which(a == uniq_ID)})))] <- 1

  mydata <- list(N1 = n, N2 = n.subj, y = patlist$nTTP,
                 X_y = X_y, W_y = W_y,
                 p1_beta_other =
                   if(dose_flag == 1){
                     c(ctrl_param$p1_beta_cycle)
                   } else{
                     c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle)
                   },
                 p2_beta_other =
                   if(dose_flag == 1){
                     matrix(ctrl_param$p2_beta_cycle, nrow = 1, ncol = 1)
                   } else {
                     diag(c(ctrl_param$p2_beta_intercept,
                            ctrl_param$p2_beta_cycle))
                   },
                 p1_beta_dose = ctrl_param$p1_beta_dose,
                 p2_beta_dose = ctrl_param$p2_beta_dose,
                 p1_tau_e = 0,
                 p2_tau_e = 1000,
                 p1_tau_g = 0,
                 p2_tau_g = 1000)

  path.model <- file.path(tempdir(), "model.file.txt")

  sink(path.model)
  cat("model
      {
      beta <- c(beta_other[], beta_dose)
      for (i in 1:N1) {
      y[i] ~ dnorm(mu[i], tau_e)
      mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], gamma)
      }
      for (j in 1:N2) {
      gamma[j] ~ dnorm(0, tau_g)
      }
      beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
      beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
      tau_e ~ dunif(p1_tau_e, p2_tau_e)           ## variance of nTTP
      tau_g ~ dunif(p1_tau_g, p2_tau_g)           ## variance of random effect
      }",fill = TRUE)
  sink()

  inits.list <- list(list(beta_other = rep(0.1, ifelse(dose_flag == 1, 1, 2)),
                          beta_dose = 0.1, tau_e = 0.1, tau_g = 0.1,
                          .RNG.seed = ceiling(runif(1) * 1e+08),
                          .RNG.name = "base::Wichmann-Hill"))

  if(length(inits.list) < n.chains){
    ## if more than 1 chain, need to modify the initial list
    inits.list <- rep(inits.list, n.chains)
    for(subc in 2:n.chains){
      inits.list[[subc]]$.RNG.seed <- ceiling(runif(1) * 1e+08) # avoid duplicate seed
    }
  }

  inits.list <- modifyList(inits.list, inits.list.set)

  jagsobj <- jags.model(path.model,
                               data = mydata,
                               n.chains = n.chains,
                               quiet = TRUE,
                               inits = inits.list)

  update(jagsobj, n.iter = burn.in, progress.bar = "none")
  post.samples <- jags.samples(jagsobj, retrieve_param,
                                      n.iter = n.iters, thin = 2,
                                      progress.bar = "none")

  return(post.samples)
}

phase1stage2 <- function(patlist, ctrl_param, n.iters = 5000, burn.in = 5000,
                         retrieve_param = c("beta_dose", "beta_other", "gamma"),
                         inits.list.set = list(), n.chains = 1, dose_flag = 0){

  ## model fitting for stage 2 (under construction) ##
  ## patlist:        the current patients record
  ## ctrl_param:     the prior setting
  ## n.iters:        number of iterations for the MCMC chain after burn in
  ## n.burn:         number of burn in
  ## retrive_param:  the parameters we want a posterior samples for the further study
  ## inits.list.mod: modification of starting point of parameters, currently use default..
  ## n.chains:       number of MCMC chains in the fitting, default is 1
  ## dose_flag:      when dose_flag equals 1, we need to degenerate design matrix
  ## return:         return the posterior samples of the MCMC chain
  #' @importFrom rjags jags.model jags.samples
  #' @importFrom stats runif update
  #' @importFrom utils modifyList

  uniq_ID <- unique(patlist$PatID)
  n.subj <- length(uniq_ID)                       ## number of patients
  n <- length(patlist$PatID)                      ## number of records
  eff.ind <- which(patlist$efficacy >= 0)         ## index of where we have efficacy record
  n.eff <- length(eff.ind)                        ## number of efficacy records
  keepeff.ind <- sapply(1:n.eff, function(i){
    which(uniq_ID == patlist$PatID[eff.ind[i]])}) ## records the location of gamma for the corresponding efficacy record

  if(dose_flag == 1){
    X_y <- cbind(patlist$cycle, patlist$dose)     ## design matrix
  } else {
    X_y <- cbind(1, patlist$cycle, patlist$dose)  ## design matrix
  }
  W_y <- matrix(0, n, n.subj)
  W_y[cbind(1:n, as.numeric(sapply(patlist$PatID, function(a) {
    which(a == uniq_ID)})))] <- 1

  y_e <- patlist$efficacy[eff.ind]                ## observed efficacy records
  X_e <- cbind(1, patlist$dose[eff.ind], patlist$dose[eff.ind]^2,
               patlist$cycle[eff.ind])
  ## design matrix of efficacy

  mydata <- list(N1 = n, N2 = n.subj, N3 = n.eff,
                 y = patlist$nTTP, X_y = X_y, W_y = W_y,
                 y_e = y_e, X_e = X_e,
                 p1_beta_other =
                   if(dose_flag == 1){
                     c(ctrl_param$p1_beta_cycle)
                   } else{
                     c(ctrl_param$p1_beta_intercept, ctrl_param$p1_beta_cycle)
                   },
                 p2_beta_other =
                   if(dose_flag == 1){
                     matrix(ctrl_param$p2_beta_cycle, nrow = 1, ncol = 1)
                   } else {
                     diag(c(ctrl_param$p2_beta_intercept,
                            ctrl_param$p2_beta_cycle))
                   },
                 p1_beta_dose = ctrl_param$p1_beta_dose,
                 p2_beta_dose = ctrl_param$p2_beta_dose,
                 p1_alpha = ctrl_param$p1_alpha,
                 p2_alpha = ctrl_param$p2_alpha,
                 p1_tau_e = 0,
                 p2_tau_e = 1000,
                 p1_tau_g = 0,
                 p2_tau_g = 1000,
                 p1_tau_f = 0,
                 p2_tau_f = 1000,
                 p1_rho = ctrl_param$p1_rho,
                 p2_rho = ctrl_param$p2_rho,
                 keepeff.ind = keepeff.ind)

  path.model <- file.path(tempdir(), "model.file.txt")
  sink(path.model)
  cat("model
      {
      beta <- c(beta_other[], beta_dose)
      for (i in 1:N1) {
      y[i] ~ dnorm(mu[i], tau_e)
      mu[i] <- inprod(X_y[i, ], beta) + inprod(W_y[i, ], gamma)
      }
      for (j in 1:N2) {
      gamma[j] ~ dnorm(0, tau_g)
      }
      for (k in 1:N3) {
      y_e[k] ~ dnorm(mu_e[k], tau_f)
      mu_e[k] <- inprod(X_e[k, ], alpha) + rho * gamma[keepeff.ind[k]]
      }

      beta_other ~ dmnorm(p1_beta_other[], p2_beta_other[, ])
      beta_dose ~ dunif(p1_beta_dose, p2_beta_dose)
      alpha ~ dmnorm(p1_alpha[], p2_alpha[, ])
      rho ~ dnorm(p1_rho, p2_rho)            ## covariance of nTTP and efficacy
      tau_e ~ dunif(p1_tau_e, p2_tau_e)      ## variance of nTTP
      tau_g ~ dunif(p1_tau_g, p2_tau_g)      ## variance of gamma
      tau_f ~ dunif(p1_tau_f, p2_tau_f)      ## variance of efficacy
      }",fill = TRUE)
  sink()
  #R2WinBUGS::write.model(model.file, path.model)
  inits.list <- list(list(beta_other = c(0.1, 0.1), beta_dose = 0.1,
                          alpha = c(0.1, 0.1, 0.1, 0.1), rho = 0.1, tau_e = 0.1,
                          tau_g = 0.1, tau_f = 0.1,
                          .RNG.seed = ceiling(runif(1) * 1e+08),
                          .RNG.name = "base::Wichmann-Hill"))

  if(length(inits.list) < n.chains){
    ## if more than 1 chain, need to modify the initial list
    inits.list <- rep(inits.list, n.chains)
    for(subc in 2:n.chains){
      inits.list[[subc]]$.RNG.seed <- ceiling(runif(1) * 1e+08) # avoid duplicate seed
    }
  }

  inits.list <- modifyList(inits.list, inits.list.set)
  jagsobj <- jags.model(path.model,
                               data = mydata,
                               n.chains = n.chains,
                               quiet = TRUE,
                               inits = inits.list)

  update(jagsobj, n.iter = burn.in, progress.bar = "none")
  post.samples <- jags.samples(jagsobj, retrieve_param,
                                      n.iter = n.iters, thin = 2,
                                      progress.bar = "none")
  return(post.samples)
}



################################################################################
####################### Function in dose recommendation  #######################
################################################################################
gen.prob.dos.cyc <- function(sim.betas, doses, cycles, thrd){

  ## Calculate the posterior probability of having safe toxicity for given
  ## doses and cycles. The safe toxicity threshold is given by thrd
  ## sim.bebas:  dim:3*N with 1st row intercept, 2nd row dose, 3rd row cycle

  M <- as.matrix(
    sapply(doses, function(d){
      sapply(cycles, function(c){
        mean(apply(sim.betas, 2, function(b){
          as.numeric(b[1] +  d * b[2] + c * b[3] <= thrd)
        }))
      })
    }))

  if (dim(M)[1] > 1 & dim(M)[2] > 1){
    rownames(M) <- paste0("C", cycles)
    colnames(M) <- paste0("D", doses)
  }
  return(M)
}


stg1.safe.dos <- function(post_samples, doses, cycles, thrd1, thrd2,
                          p_tox1, p_tox2, ICD.flag){

  ## decide the safe dose after model fittng ##

  sim.betas <- as.matrix(rbind(post_samples$beta_other[1, , 1],
                               post_samples$beta_dose[, , 1],
                               post_samples$beta_other[2, , 1]))

  prob.dose.cycle1 <- gen.prob.dos.cyc(sim.betas = sim.betas,
                                       doses = doses, cycles = 1,
                                       thrd = thrd1)

  if(ICD.flag == T){
    # no longer need cycle 2-6 crition if ICD.flag == T
    allow.doses <- which(prob.dose.cycle1[, 1] >= p_tox1)
  } else {
    # check the whether the dose is safe across subsequent cycle if ICD.flag = F
    prob.dose.cycle2 <- gen.prob.dos.cyc(sim.betas = sim.betas,
                                         doses = doses, cycles = cycles[-1],
                                         thrd = thrd2)
    allow.doses <- which(prob.dose.cycle1[, 1] >= p_tox1 &
                           colMeans(prob.dose.cycle2) >= p_tox2)
  }
  return(allow.doses)
}

stg1.dos.rec.i <- function(post_samples, uniq_ID, patID_act,
                           cycle_act, rec_dose_act, Max_tested_doseA,
                           doses, cycles, c1, p1, DLT.drop.flag = F, y.dlt = NULL,
                           testedD = T){

  ## Calculate the posterior probability of having safe toxicity for each patients
  ## for each doses and cycles. The safe toxicity threshold is given by thrd.
  ## And recommend dose for the activtive patients for the next cycle

  ## post_samples:  posteior samples from stage 1 model fitting
  ## uniq_ID:       ID for all enrolled patients
  ## patID_act:     current active patients' ID
  ## cycle_act:     current cycle for active patients
  ## rec_dose_act:  current recommend dose for current cycle
  ## Max_tested_doseA:      Maximum tested dose level in cycle1
  ## doses:         dose levels in the study
  ## cycles:        cycles in the treatment
  ## DLT.drop.flag: whether drop off when DLT observed
  ## y.dlt:         the DLT results, required when DLT.drop.flag = T
  ## c1:            target toxicity
  ## p1:            ICD_thrd
  ## testedD:       whether recommend dose for new cohort among cycle 1 tested dose level

  sim.post <- as.matrix(rbind(post_samples$beta_other[1, , 1],
                              post_samples$beta_dose[, , 1],
                              post_samples$beta_other[2, , 1],
                              post_samples$gamma[, , 1]))

  # find the activie patients need prediction for next cycle
  cycle_nxt <- cycle_act + 1
  if(DLT.drop.flag == T){
    nxt.index <- which(cycle_nxt <= max(cycles) & y.dlt == 0)
  } else {
    nxt.index <- which(cycle_nxt <= max(cycles))
  }
  cycle_nxt <- cycle_nxt[nxt.index]
  patID_nxt <- patID_act[nxt.index]

  if(length(cycle_nxt) == 0){
    return(list(patID_nxt = NULL, cycle_nxt = NULL, rec_dose_nxt = NULL,
                ICD_nxt = NULL))}

  n.subs <- length(patID_nxt)
  Pat.index <-as.numeric(sapply(patID_nxt, function(a) {
    which(a == uniq_ID)}))

  ## M store the posterior probability
  M <-
    sapply(1:n.subs, function(i){
      sapply(doses, function(d){
        mean(apply(sim.post, 2, function(b){
          as.numeric(b[1] +  d * b[2] + cycle_nxt[i] * b[3] +
                       b[Pat.index[i] + 3] <= c1)
        }))
      })
    })

  ## find the maximum dose that can be consider as safe and use as recommend dose for next level
  rec_dose_nxt <- apply(M, 2, function(x){max(which(x >= p1))})

  if(testedD == T){
    ## don't escalcate dose to untested dose level ##
    rec_dose_nxt[which(rec_dose_nxt > Max_tested_doseA)] <- Max_tested_doseA
  }
  ## don't skip dose level ##
  skip.index <- which((rec_dose_nxt - rec_dose_act[nxt.index]) > 1)
  rec_dose_nxt[skip.index] <- rec_dose_act[nxt.index[skip.index]] + 1

  ## remove patients who has no safe dose level for next cycle
  hazard.index <- which(rec_dose_nxt < -1)
  if(length(hazard.index) > 0){
    cycle_nxt <- cycle_nxt[-hazard.index]
    patID_nxt <- patID_nxt[-hazard.index]
    n.subs <- length(patID_nxt)
    rec_dose_nxt <- rec_dose_nxt[-hazard.index]
  }

  if(length(cycle_nxt) == 0){
    return(list(patID_nxt = NULL, cycle_nxt = NULL, rec_dose_nxt = NULL,
                ICD_nxt = NULL))}

  ## check the ICD information, just used in developing packages ##
  ICD_nxt <- paste0("ICD: D", rec_dose_nxt, ": ",
                    format(round(M[cbind(rec_dose_nxt, 1:n.subs)], 3),
                           nsmall= 3))

  return(list(patID_nxt = patID_nxt, cycle_nxt = cycle_nxt,
              rec_dose_nxt = rec_dose_nxt,
              ICD_nxt = ICD_nxt))
}

stg2.eff.rec.i <- function(post_samples, cycle_nxt, rec_dose_nxt, proxy.thrd){

  ## after finding the safe dose for each active patients for next cycle,
  ## calculate the expected population efficacy across all the safe doses and cycle
  ## then recommend dose for the activtive patients for the next cycle
  ## need to run stg1.dos.rec.i first
  ## proxy.thrd: the buffer of choosing IED

  if(length(cycle_nxt) == 0){return(list(rec_dose_nxt = NULL, IED_nxt = NULL))}

  sim.alphas <- as.matrix(rbind(post_samples$alpha[, , 1]))

  # calculate the dictionary
  Pop.EFF <- sapply(1:max(rec_dose_nxt), function(d){
    sapply(1:max(cycle_nxt), function(c){
      mean(apply(sim.alphas, 2, function(o){
        as.numeric(o[1] + o[2] * d + o[3] * d^2 + o[4] * c)
      }))
    })
  })

  rec_dose_nxt <- apply(cbind(cycle_nxt, rec_dose_nxt), 1,
                        function(a){
                          which(Pop.EFF[a[1], 1:a[2]] >=
                                  (max(Pop.EFF[a[1], 1:a[2]]) - proxy.thrd))[1]})
  # we tend to choose lower dose level when the posterior predictive estimate of efficacy
  # is no less than Max.effecy - proxy.thrd

  ## check the ICD information, just used in developing packages ##
  IED_nxt <- paste0("IED: D", rec_dose_nxt, ": ",
                    format(round(Pop.EFF[cbind(cycle_nxt, rec_dose_nxt)], 3),
                           nsmall= 3))

  return(list(rec_dose_nxt = rec_dose_nxt, IED_nxt = IED_nxt))
}

stg2.eff.rec.i.all.cycle <- function(post_samples, allow.doses, cycles,
                                     proxy.thrd){

  ## after finding the safe dose of the trial,
  ## calculate efficacious dose among the safe dose for all cycle
  ## allow.doses: the safe dose of the trial
  ## cycles: the treatment cycles of the trial
  ## proxy.thrd: the buffer of choosing efficacious dose
  ## return:
  ## recom.MED.cycle: recommended MED for all cycles

  sim.alphas <- as.matrix(rbind(post_samples$alpha[, , 1]))

  ## calculate posterior estimation of efficacy for all allowable dose all cycles
  PEFF <- sapply(allow.doses, function(d) {
    sapply(cycles, function(c){
      mean(apply(sim.alphas, 2, function(o) {
        as.numeric(o[1] + o[2] * d + o[3] * d^2 + o[4] * c)    # estimate efficacy for dose a cycle 1
      }))
    })
  })

  recom.MED.cycle <- apply(PEFF, 1, function(o){
    cycles[which(o >= max(o) - proxy.thrd)[1]]
  })

  return(list(recom.MED.cycle = recom.MED.cycle))
}


stg1.dos.rec <- function(post_samples, doses, tox.target, Max_tested_doseA){

  ## dose recommendation for next cohort for stage 1 ##
  ## doses
  ## tox.target = 0.28
  ## Max_tested_doseA: the maximum previously tried dose for cycle1

  sim.betas <- as.matrix(rbind(post_samples$beta_other[1, , 1],
                               post_samples$beta_dose[, , 1],
                               post_samples$beta_other[2, , 1]))

  loss.doses <- sapply(doses, function(d) {
    mean(apply(sim.betas, 2, function(b) {
      abs(b[1] + b[2] * d + b[3] - tox.target)
    }))
  })

  nxtdose <- doses[which.min(loss.doses)]

  # No skipping of dose levels not previously tried is allowed, so if there happens to be skipping of dose
  #levels from minimizing the Bayesian risk, we will de-escalate to the next safer dose that does not
  # skip any dose levels not previously tried.
  if (as.numeric(nxtdose) > Max_tested_doseA + 1) {
    doseA <- Max_tested_doseA + 1
  } else {
    doseA <- as.numeric(nxtdose)
  }

  return(doseA)
}

stg3.dos.rec <- function(post_samples, doses, tox.target, Max_tested_doseA){

  ## dose recommendation for next cohort for stage 1 ##
  ## doses
  ## tox.target = 0.28
  ## Max_tested_doseA: the maximum previously tried dose for cycle1

  sim.betas <- as.matrix(rbind(post_samples$beta_other[1, , 1],
                               post_samples$beta_dose[, , 1],
                               post_samples$beta_other[2, , 1]))

  loss.doses <- sapply(doses, function(d) {
    mean(apply(sim.betas, 2, function(b) {
      abs(b[1] + b[2] * d + b[3] - tox.target)
    }))
  })

  nxtdose <- doses[which.min(loss.doses)]

  # Don't recommend dose levels not tested in cycle 1 in stage 3
  if (as.numeric(nxtdose) > Max_tested_doseA ) {
    doseA <- Max_tested_doseA
  } else {
    doseA <- as.numeric(nxtdose)
  }

  return(doseA)
}

stg2.eff.rec <- function(post_samples, allow.doses, Max_tested_doseA,
                         proxy.thrd){

  ## dose recommendation for next cohort for stage 2(with efficacy data) ##
  ## doses
  ## tox.target = 0.28
  ## Max_tested_doseA: the maximum previously tried dose for cycle1
  ## proxy.thrd: the buffer of efficacy

  sim.alphas <- as.matrix(rbind(post_samples$alpha[, , 1]))

  ## calculate posterior estimation of efficacy for all allowable dose
  RAND.EFF <- sapply(allow.doses, function(a) {
    mean(apply(sim.alphas, 2, function(o) {
      as.numeric(o[1] + o[2] * a + o[3] * a^2 + o[4])    # estimate efficacy for dose a cycle 1
    }))
  })

  ## randomly sample next recommended dose
  #  RAND.EFF <- exp(RAND.EFF) / sum(exp(RAND.EFF))
  #  nxtdose <- sample(allow.doses, 1, prob = RAND.EFF)

  ## next dose is the lowest dose (among allow dose set) that is efficacious
  nxtdose <- allow.doses[which(RAND.EFF >= max(RAND.EFF) - proxy.thrd)[1]]

  # No skipping of dose levels not previously tried is allowed, so if there happens to be skipping of dose
  #levels from minimizing the Bayesian risk, we will de-escalate to the next safer dose that does not
  # skip any dose levels not previously tried.
  if (nxtdose > Max_tested_doseA + 1) {
    nxtdose <- Max_tested_doseA + 1
  }
  return(nxtdose)
}


stg3.eff.rec <- function(post_samples, allow.doses, Max_tested_doseA,
                         proxy.thrd){

  ## Find the recommended dose for the end of the study (stage3)
  # proxy.thrd: The allowed proxy bandwidth for efficacy

  sim.alphas <- as.matrix(rbind(post_samples$alpha[, , 1]))

  ## calculate posterior estimation of efficacy for all allowable dose for cycle 1
  effcy.dose <- sapply(allow.doses, function(a) {
    mean(apply(sim.alphas, 2, function(o) {
      as.numeric(o[1] + o[2] * a + o[3] * a^2 + o[4])
    }))
  })

  proxy.eff.doses <- allow.doses[which(sapply(effcy.dose, function(a){
    abs(a - max(effcy.dose)) <= proxy.thrd
  }))]
  recom.doses <- min(proxy.eff.doses)

  # Don't recommend untested dose level
  if (recom.doses > Max_tested_doseA) {
    recom.doses <- Max_tested_doseA
  }

  return(recom.doses)
}


################################################################################
###########################  Function in diagonosis   ##########################
################################################################################

pp.nTTP <- function(post_samples, doses, cycles, target){

  ## calculate the posterior probability of nttp < target for all cycles and doses
  # post_samples: posterior samples of stage3 model fitting
  # doses: all dose levels
  # cycles: all cycles
  # target toxicity

  sim.betas <- as.matrix(rbind(post_samples$beta_other[1, , 1],
                               post_samples$beta_dose[, , 1],
                               post_samples$beta_other[2, , 1]))

  pp.nTTPM <- sapply(doses, function(d){
    sapply(cycles, function(c){
      mean(apply(sim.betas, 2,
                 function(b){
                   as.numeric(b[1] +  d * b[2] + c * b[3] <= target)
                 }))
    })
  })
  colnames(pp.nTTPM) <- paste0("D", doses)
  rownames(pp.nTTPM) <- paste0("C", cycles)
  return(pp.nTTPM)
}

dlt.rt.c1 <- function(list_simul, chSize){
  # calculate cycle 1 dlt rate
  sum(sapply(list_simul, function(a){
    sum(a$patlist$cycle[which(a$patlist$dlt == 1)] == 1)})) /
    (sum(sapply(list_simul, function(a){a$n.cohort})) * chSize)
}

dlt.rt.subseq <- function(list_simul, chSize){
  # calculate dlt rate on cycle > 1
  sum(sapply(list_simul, function(a){
    sum(a$patlist$cycle[which(a$patlist$dlt == 1)] > 1)})) /
    (sum(sapply(list_simul, function(a){a$n.cohort})) * chSize)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  ## Summarizes data.
  ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   groupvars: a vector containing names of columns that contain grouping variables
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  #' @importFrom plyr ddply rename

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, sd and quantiles
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     q.2.5 = as.numeric(quantile(xx[[col]], na.rm = na.rm, probs = 0.025)),
                     q.97.5 = as.numeric(quantile(xx[[col]], na.rm = na.rm, probs = 0.975))
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}





