SimPRMD <-function(seed = 1234, numTrials = 100, doses = 1:6, cycles = 1:6,
                   eff.structure = matrix(0, nrow = 6, ncol = 6),
                   eff.Sigma = diag(6), eff.sd_trans = 1.5, tox.target= 0.28,
                   p_tox1 = 0.2, p_tox2 = 0.2, trialSize = 36, chSize = 3,
                   thrd1 = 0.28, thrd2 = 0.28, proxy.thrd = 0.1,
                   tox.matrix = NULL,
                   wm = matrix(c(0, 0.5, 0.75, 1, 1.5,
                                 0, 0.5, 0.75, 1, 1.5,
                                 0, 0, 0, 0.5, 1),
                               byrow = T, ncol = 5), toxmax = 2.5,
                   toxtype = NULL, intercept.alpha = NULL,
                   coef.beta = NULL, cycle.gamma = NULL,
                   param.ctrl = list(),
                   n.iters = 10000, burn.in = 5000, thin = 2,
                   n.chains = 1, effcy.flag = T, ICD.flag = T,
                   DLT.drop.flag = T, testedD = T,  IED.flag = T,
                   ICD_thrd = 0.3) {

  #' Simulation for a Multi-Stage Phase I Dose-Finding Design
  #'
  #' A function to implement simulations for a multi-stage phase 1 dose-finding
  #' design incorporating a longitudinal continuous efficacy outcome and
  #' toxicity data from multiple treatment cycles. The available models include
  #' 1-stage model with/without individualized dose modification, 3-stage model
  #' with/without individualized dose modification, 3-stage model with
  #' individualized dose modification on stage II and 3-stage model with
  #' individualized dose modification on stage I and dose modification on stage
  #' II.
  #'
  #' @param seed The seed of R's random number generator. Default is 1234
  #' @param numTrials An integer specifying the number of simulations
  #' @param doses A vector of doses that users are going to explore. Default is
  #'   1:6, where dose 1 through dose 6 are being tested.
  #' @param cycles A vector of cycles that the treatment plans to go through.
  #'   Default is 1:6, where patients will experience up to 6 cycles of the
  #'   treatment
  #' @param eff.structure  A matrix provides the mean of the multivariate
  #'   Gaussian distribution in efficacy data generation. Specifically, the
  #'   \eqn{(i, j)}th element represents the mean value of \eqn{i}th dose level
  #'   and \eqn{j}th cycle of the Gaussian distribution for efficacy data
  #'   generation. Default is a 6 by 6 zero matrix
  #' @param eff.Sigma The covariance matrix of the multivariate Guassian
  #'   distribution in efficacy data generation. See details below.
  #' @param eff.sd_trans A positive number controls the skewness of the
  #'   distribution of the efficacy response. Default is 1.5. See details below.
  #' @param tox.target  The target toxicity of the treatment. Default is 0.28.
  #'   See details below.
  #' @param p_tox1  The probability cutoff for cycle 1 toxicity. Default is 0.2.
  #'   See details below.
  #' @param p_tox2  The probability cutoff for later cycles toxicity beyond
  #'   cycle 1. Default is 0.2. See Details below.
  #' @param trialSize The maximum sample size for trial simulation. Default is
  #'   36. Must be the multiple of cohort size, represented by chSize
  #' @param chSize    The cohort size of patients recruited. Default is 3.
  #' @param thrd1 An upper bound of toxicity for cycle 1 of the treatment.
  #'   Default is 0.28. See Details below.
  #' @param thrd2 An upper bound of toxicity for late cycles of the treatment,
  #'   beyond cycle 1. Default is 0.28. See Details below
  #' @param proxy.thrd A distance parameter to define efficacious doses. Any
  #'   dose whose predicted efficacy is within proxy.thrd away from the largest
  #'   one among the safe doses will be declared an efficacious dose.
  #' @param tox.matrix Optional. A four-dimension array specifying the
  #'   probabilities of the occurrences of certain grades for certain types of
  #'   toxicities, at each dose level and cycle under consideration. Dimension 1
  #'   refers to doses; dimension 2 corresponds to cycles of the treatment;
  #'   dimension 3 regards the types of toxicities while dimension 4 relates to
  #'   grades. If null, which is default choice, the arguments toxtype,
  #'   intercept.alpha, coef.beta, cycle.gamma must be provided to simulate this
  #'   array.
  #' @param wm Clinical weight matrix, where toxicity types define the rows
  #'   while the toxicity grades define the columns. Usually solicited from
  #'   physicians.
  #' @param toxmax The normalization constant used in computing nTTP score. For
  #'   details, see Ezzalfani et al(2013).
  #' @param toxtype Only specified when tox.matrix is null. This argument, a
  #'   character vector, specifies toxicity types considered in the trial.
  #' @param intercept.alpha Only specified when tox.matrix is null. A four
  #'   element numeric vector specifying the intercepts for the cumulative
  #'   probabilities of the occurrences of grades 0-4 of toxicities in
  #'   proportional odds model. See Details below.
  #' @param coef.beta Only specified when tox.matrix is null. A n numeric vector
  #'   specifying the slope for dose in proportional odds model for n types of
  #'   toxicities. See Details below
  #' @param cycle.gamma Only specified when tox.matrix is null. A scalar
  #'   controlling the cycle effect in simulation in proportional odds model.
  #'   See Details below
  #' @param param.ctrl A list specifying the prior distribution for the
  #'   parameters. \describe{\item{p1_beta_intercept}{the prior mean of
  #'   intercept of toxicity model assuming a normal prior}
  #'
  #'   \item{p2_beta_intercept}{the precision (inverse of variance) of intercept
  #'   of toxicity model assuming a normal prior}
  #'
  #'   \item{p1_beta_cycle}{the prior mean of cycle effect of toxicity model
  #'   assuming a normal prior}
  #'
  #'   \item{p2_beta_cycle}{the precision (inverse of variance) of cycle effect
  #'   of toxicity model assuming a normal prior}
  #'
  #'   \item{p1_beta_dose}{the prior minimum of dose effect of toxicity model
  #'   assuming a uniform prior}
  #'
  #'   \item{p2_beta_dose}{the prior maximum of dose effect of toxicity model
  #'   assuming a uniform prior}
  #'
  #'   \item{p1_alpha}{the prior mean vector of the parameters from efficacy
  #'   model assuming a multivariate normal prior}
  #'
  #'   \item{p2_alpha}{the prior precision matrix (inverse of covariance matrix)
  #'   of the parameters from efficacy model assuming a multivariate normal
  #'   prior}
  #'
  #'   \item{p1_gamma0}{the prior mean of association parameter \eqn{\gamma}
  #'   (See Du et al(2017)) of two submodels of the joint model assuming a
  #'   normal prior}
  #'
  #'   \item{p2_gamma0}{the prior precision (inverse of variance) of association
  #'   parameter \eqn{\gamma} of two submodels of the joint model assuming a
  #'   normal prior. } Default is non-informative priors. }
  #' @param n.iters Total number of MCMC simulations. Default is 10,000.
  #' @param burn.in Number of burn=ins in the MCMC simulation. Default is 5,000.
  #' @param thin Thinning parameter. Default is 2.
  #' @param n.chains  No. of MCMC chains in Bayesian model fitting. Default is 1
  #' @param DLT.drop.flag Whether the patients should suspend the treatment when
  #'   observing DLT. Default is TRUE
  #' @param effcy.flag Whether we include efficacy response in modeling or not?
  #' @param ICD.flag Whether we allow dose changing for cycle > 1 in stage 1
  #'   model or not? Default is TRUE. See details below
  #' @param testedD Default is TRUE. Whether we only allow ICD or IED among
  #'   cycle 1 tested dose level
  #' @param IED.flag Default is TRUE. Whether we allow dose changing for cycle >
  #'   1 in stage 2 model or not?
  #' @param ICD_thrd The cut-off point of the posterior toxicity probability in
  #'   defining ICD. Default is 0.3. See details below.
  #'
  #' @details The user can simulation efficacy response with different
  #'   dose-efficacy and cycle-efficacy pattern using argument
  #'   \code{eff.structure}, \code{eff.Sigma} and \code{eff.sd_trans}. The
  #'   sampling process of efficacy response start from generating sample \eqn{z
  #'   = {z1, \ldots, zd} } from multivariate Gaussian distribution \deqn{z ~
  #'   MVN(\mu, V)}, where \eqn{\mu} and \eqn{V} are specified by
  #'   \code{eff.structure} and \code{eff.Sigma}, respectively. Define
  #'   \eqn{\phi} be the density of \eqn{N(0, \sigma^2)} with CDF \eqn{\Phi},
  #'   and \eqn{\sigma^2} is set by \code{eff.sd_trans}. Then the efficacy
  #'   response is calculated by taking the CDF of \eqn{z}: \deqn{x={x1, \ldots,
  #'   xd} = \Phi(z) = { \Phi(z1), \ldots, \Phi(zd)}} is the generated efficacy
  #'   response. Notice here the variance parameter \eqn{\sigma^2_{trans}}
  #'   controls the variance of the generated efficacy.
  #'
  #' @return \item{senerio_sum}{contains \code{mnTTP.M} the matrix of mean nTTP
  #'   for each dose and cycle and \code{pDLT.M} matrix of probability of
  #'   observing DLT for each dose and cycle} \item{eff_sum}{When
  #'   \code{effcy.flag == TRUE}, contains \code{eff.M} the mean efficacy for
  #'   each dose and cycle and \code{err.cor.ls} A list with a length of dose
  #'   levels numbers recording the marginal correlation matrix across cycles of
  #'   efficacy data for each dose level} \item{list_simul}{A list of length
  #'   numTrials. Each element includes \code{patlist} which records all the
  #'   treatment and outcome information; \code{dose_aloca} which shows the
  #'   cycle 1 dose allocation; \code{doseA} which saves the recommended dose
  #'   level for cycle 1 at the end of the phase I simulation, equals "early
  #'   break" if the trial was stop before finishing the trial; \code{n.cohort}
  #'   indicates the last cohort in the trial; \code{pp.nTTPM} gives the
  #'   posterior probability of nTTP less than target toxicity \code{tox.target}
  #'   for all dose level any cycles and \code{message} saves the message of
  #'   each trial.} \item{chSize}{The input argument \code{chSize}}
  #'   \item{sim.time}{Time cost in simulation} \item{doses}{The input argument
  #'   \code{doese}} \item{cycles}{The input argument \code{cycles}}
  #'   \item{effcy.flag}{The input argument \code{effcy.flag}}
  #'   \item{proxy.thrd}{The input argument \code{proxy.thrd}}
  #'   \item{DLT.drop.flag}{The input argument \code{DLT.drop.flag}}
  #'
  #' @details The user can simulate longitudinal efficacy response with
  #'   different dose-efficacy and cycle-efficacy pattern using argument
  #'   \code{eff.structure}, \code{eff.Sigma} and \code{eff.sd_trans}. The
  #'   sampling process of efficacy response starts from generating \eqn{z =
  #'   {z1, \ldots, zd} } from multivariate Gaussian distribution \deqn{z ~
  #'   MVN(\mu, V)}, where \eqn{\mu} and \eqn{V} are specified by
  #'   \code{eff.structure} and \code{eff.Sigma}, respectively. Define
  #'   \eqn{\phi} be the density of \eqn{N(0, \sigma^2)} with CDF \eqn{\Phi},
  #'   where \eqn{\sigma^2} is set by \code{eff.sd_trans}. Then the efficacy
  #'   measure is generated by taking the CDF of \eqn{z}: \deqn{x={x1, \ldots,
  #'   xd} = \Phi(z) = { \Phi(z1), \ldots, \Phi(zd)}}. Notice here the variance
  #'   parameter \eqn{\sigma^2_{trans}} controls the variance of the generated
  #'   efficacy.
  #'
  #'   \code{p_tox1}, \code{p_tox2}, \code{thrd1} and \code{thrd2} are used to
  #'   define allowable (safe) doses the probability conditions for cycle 1:
  #'   \deqn{P(nTTP1 < thrd1) > p_tox1} and for cycle > 1: \deqn{p(nTTP2 <
  #'   thrd2) > p_tox2} , where \eqn{nTTP1} and \eqn{nTTP2} denote the posterior
  #'   estimate of nTTP for cycle 1 and the average of cycle > 1. When we
  #'   implement model with individualized dose modification, we only check the
  #'   condition for cycle 1 for defining allowable (safe) doses.
  #'
  #'   \code{ICD_thrd} are used to find ICD. ICD is defined as the maximum dose
  #'   which satisfy the condition \deqn{P(nTTPi <  target.tox) > ICD_thrd} ,
  #'   where \eqn{nTTPi} is the individualized posterior predicted nTTP score.
  #'   The individualized dose modification for next cycle will not escalate
  #'   more than 1 dose from the current dose.
  #'
  #'
  #' @examples
  #' data("prob")      # load prob.RData from package phaseI, Details see "?prob"
  #' data("eff")       # load eff.RData from package phaseI. Details see "?eff"
  #'
  #' eff.structure = eff$Dose_Cycle_Meff[2, 2, , ]
  #' eff.Sigma = eff$Sigma
  #' eff.sd_trans = eff$sd_trans
  #'
  #' wm <- matrix(c(0, 0.5, 0.75, 1, 1.5,
  #'                0, 0.5, 0.75, 1, 1.5,
  #'                0, 0, 0, 0.5, 1),
  #'              byrow = TRUE, ncol
  #'               = 5)                          # weighted matrix for toxicity matrix
  #'                                             # nrow = No.of type; ncol = No. of grade
  #' toxmax <- 2.5
  #' tox.matrix <- prob["MTD4", "flat", , , , ]
  #'
  #'
  #' #------- a flat dose-toxicity, dose-efficacy, cycle-efficacy pattern------#
  #' \donttest{
  #' simul1 <- SimPRMD(numTrials = 1, tox.matrix = tox.matrix,
  #'                   eff.structure = eff.structure, eff.Sigma = eff.Sigma,
  #'                   eff.sd_trans = eff.sd_trans, wm = wm, toxmax = toxmax,
  #'                   trialSize = 36)
  #' }
  #' #------- a flat dose-toxicity pattern model ------#
  #' \donttest{
  #' simul2 <- SimPRMD(numTrials = 1, toxtype = c("H", "L", "M"),
  #'                   intercept.alpha = c(1.9, 2.3, 2.6, 3.1),
  #'                   coef.beta = c(-0.3, -0.2, -0.25),
  #'                   cycle.gamma = 0, tox.target = 0.23,
  #'                   thrd1 = 0.23, thrd2 = 0.23, p_tox1 = 0.2, p_tox2 = 0.2,
  #'                   ICD.flag = FALSE, IED.flag = FALSE, effcy.flag = TRUE)
  #'
  #' summary(simul2)
  #' plot(simul2)
  #' }
  #'
  #' @import coda
  #' @importFrom arrayhelpers vec2array
  #' @importFrom phase1RMD GenToxProb
  #' @importFrom coda gelman.diag
  #' @importFrom utils setTxtProgressBar
  #' @importFrom utils txtProgressBar
  #' @export

  #######################################################################
  # input argments checking (modified later)
  #######################################################################

  if (nrow(eff.structure) != length(doses) |
      ncol(eff.structure) != length(cycles)){
    stop("Check if you have specified the efficacy mean structure for the correct number of doses. \n")
  }

  if(is.null(wm)){
    stop("The clinical weight matrix for toxicities is not specified \n")
  }

  ## test whether trialSize/chSize is integer
  if(trialSize %% chSize != 0){
    stop("trialSize has to be the multiple of cohort size\n")
  }

  if(effcy.flag == T & is.null(proxy.thrd)){
    stop("proxy.thrd is required when effcy.flag == TRUE\n")
  }
  #######################################################################
  # Model fitting report
  #######################################################################

  if(effcy.flag == F){
    if(ICD.flag == T){
      cat("1-stage model with individualized dose modification (ICD on) \n")
    } else{
      cat("1-stage model (ICD off) \n")
    }
  } else {
    if(IED.flag == F){
      if(ICD.flag == T){
        cat("3-stage model with individualized dose modification in stage I and ")
        cat("dose \nmodification in stage II (ICD on, IED off) \n")
      }else{
        cat("3-stage model (ICD off, IED off) \n")
      }
    } else {
      if(ICD.flag == T){
        cat("3-stage model with individualized dose modification (ICD on, IED on) \n")
      } else {
        cat("3-stage model with individualized dose modification only in stage II (ICD off in stage I, IED on) \n")
      }
    }
  }
  if(DLT.drop.flag == T){
    cat("Patients will not continue treatment when having DLT \n")
  }

  ########################################################################
  # precalculate data
  ########################################################################

  MaxCycle <- length(cycles)
  Numdose <- length(doses)
  listdo <- paste0("D", doses)
  nTTP.all <- nTTP.array(wm, toxmax)   ## generate array of nTTP for later usage
  inits.list.set <- list()
  if(is.null(tox.matrix)){
    ## Generate the tox prob matrix if the tox.matrix is not given
    ## Need to be modified
    tox.matrix = GenToxProb(toxtype, intercept.alpha, coef.beta,
                            cycle.gamma)
  }

  # obtain the mnTTP.M and pDLT.M for the senerio
  senerio_sum <- nTTP_summary(tox.matrix, nTTP.all, wm)

  if(effcy.flag == T){
    # obtain the eff.M for the senerio
    eff_sum <- eff_summary(eff.structure = eff.structure, eff.Sigma = eff.Sigma,
                           eff.sd_trans = eff.sd_trans, seed = seed,
                           plot.flag = F)
  }else{
    eff_sum <- NULL
  }

  # default prior settings for all model fitting #
  ctrl_param <- list(p1_beta_intercept = 0, p1_beta_cycle = 0,
                     p2_beta_intercept = 0.001, p2_beta_cycle = 0.001,
                     p1_beta_dose = 0, p2_beta_dose = 1000,
                     p1_alpha = c(0, 0, 0, 0), p2_alpha = diag(rep(0.001, 4)),
                     p1_rho = 0, p2_rho = 0.001)

  ## modifty setting if specific setting is given
  ctrl_param <- modifyList(ctrl_param, param.ctrl)


  ########################################################################
  # simulate data
  ########################################################################

  t <- proc.time()
  set.seed(seed)
  list_simul <- list()
  seed_rand <- ceiling(runif(numTrials) * 1000000000)
  cat("Simulation process:\n")
  pb <- txtProgressBar(min = 0, max = numTrials, style = 3)
  for(i in 1:numTrials) {
    set.seed(seed_rand[i])
    ## pat_list: records of simuated patients
    patlist <- list(PatID = NULL, dose = NULL, cycle = NULL,
                    nTTP = NULL, dlt = NULL, efficacy = NULL, effz= NULL)

    patID_act<- NULL                          ## patID_list:   the ID of the active patients in the study
    rec_dose_act <- NULL                      ## rec_dose_act: recommend dose for the active patient in the study
    cycle_act <- NULL                         ## cycle_act:    the current cycle for the active patient
    dose_aloca <- rep(0, Numdose)             ## dose_aloca:   dose allocation record for cycle 1

    ## for the begining of the phase I trill
    doseA <- min(doses)                       ## doseA:        the dose assigned for the 1st cohort
    rec_doseA <- NULL                         ## rec_doseA:    record all the assigned dose for the next cohort
    Max.cohort <- trialSize / chSize          ## maximum number of cohort
    end.label <- FALSE
    break.label <- FALSE                      ## whether the simulation is breaked before stage 3
    message <- NULL                           ## record the message of each simulation, default is NULL.
    n.cohort <- 1
    while(n.cohort < (Max.cohort + 1) | length(cycle_act) > 0){
      ## end the simulation when no new cohort and no active patients

      ############################
      # Generate simulation data #
      ############################

      if(n.cohort <= Max.cohort ){

        # generate ID for new recruited patients
        patID_act <- c(patID_act, paste0("cohort", n.cohort, "subject", 1:chSize))

        # generate dose for all active patients
        rec_dose_act <- c(rec_dose_act, rep(doseA, chSize))

        # generate cycle for all active patients
        cycle_act <- c(cycle_act, rep(1, chSize))

        # record dose allocation and new assigned dose for new cohort
        dose_aloca[doseA] <- dose_aloca[doseA] + chSize
        rec_doseA <- c(rec_doseA, doseA)

      }

      # generate the outcome for active patients and save the record
      outcome <- apply(cbind(rec_dose_act, cycle_act, patID_act), 1,
                       gen_nTTP_dlt, tox.matrix, wm, nTTP.all,
                       eff.structure, eff.Sigma, eff.sd_trans, patlist)


      # updata patient, dose and cycle in the dataset
      n.point <- length(patlist$PatID)
      patlist$PatID[(n.point + 1):(n.point + length(patID_act))] <- patID_act
      patlist$dose[(n.point + 1):(n.point + length(patID_act))] <- rec_dose_act
      patlist$cycle[(n.point + 1):(n.point + length(patID_act))] <- cycle_act

      patlist$nTTP <- c(patlist$nTTP, outcome["y.nTTP", ])
      patlist$dlt <- c(patlist$dlt, outcome["y.dlt", ])
      patlist$efficacy <- c(patlist$efficacy, outcome["y.eff", ])
      patlist$effz <- c(patlist$effz, outcome["y.effz", ])

      ######################################
      # model fitting & dose recommendation#
      ######################################

      ## rec_dose_act:              recommended dose for the active patient for next cycle
      ## doseA:                     recommended dose for the new cohort

      if (n.cohort == 1) {

        ## the first cohort, use 3 + 3 design

        dlt <- outcome["y.dlt", ]     ## new generated dlt data
        if (sum(dlt) == 0) {
          doseA <- doseA + 1          ## no dlt observed, escalate dose to 2
          dose_flag <- 0              ## whether the current dose observe dlt
        } else if (sum(dlt <= 2)){
          doseA <- doseA              ## observe 1 or 2 dlts, same dose for cohort 2
          dose_flag <- 1
        } else{
          message <- paste(message, "\n early stop in the first cohort, 3 dlt observed")
          break.label = T
          break
        }

        if(DLT.drop.flag == T){
          # DLT drop #
          nxt.index <- which(outcome["y.dlt", ] == 0)
          patID_act <- patID_act[nxt.index]
          rec_dose_act <- rec_dose_act[nxt.index]  ## don't change dose in 3+3 design
          cycle_act <- cycle_act[nxt.index] + 1    ## update cycle
        } else{
          cycle_act <- cycle_act + 1    ## update cycle
        }
      } else if((n.cohort == 2) & (dose_flag == 1)){

        ## the second cohort, use 3 + 3 design

        dlt <- patlist$dlt[which(patlist$cycle == 1)]
        ## 6 records of 2 cohort in the first cycle
        if (sum(dlt) == 1) {
          doseA <- doseA + 1          ## 1 dlt out of 6 records, escalate to dose 2
          dose_flag <- 0              ## change flag to 0
        } else if(sum(dlt) == 2){
          doseA <- doseA              ## 2 dlts out of 6 records, same dose
          dose_flag <- 1
        } else {
          message <- paste(message, "\n early stop, ",
                           "more than two 2 dlts out of 6 records",
                           "in the first 2 cohort")
          break.label = T
          break
        }

        if(DLT.drop.flag == T){
          # DLT drop #
          nxt.index <- which(outcome["y.dlt", ] == 0)
          patID_act <- patID_act[nxt.index]
          rec_dose_act <- rec_dose_act[nxt.index]  ## don't change dose in 3+3 design
          cycle_act <- cycle_act[nxt.index] + 1    ## update cycle
        } else{
          cycle_act <- cycle_act + 1    ## update cycle
        }
      } else if(ifelse(effcy.flag == T,         # If efficacy is true
                       n.cohort < (Max.cohort / 2),
                       n.cohort <= (Max.cohort + 1))){

        ### stage 1 ###
        ## model fitting ##
        ## model fitting for n.cohorts = 5 (<6) in the example ##
        if(n.cohort < Max.cohort | ICD.flag == T){
          retrieve_param <- c("beta_dose", "beta_other", "gamma")
          post_samples <- phase1stage1(patlist = patlist,
                                       ctrl_param = ctrl_param,
                                       n.iters = n.iters - burn.in,
                                       burn.in = burn.in,
                                       retrieve_param = retrieve_param,
                                       dose_flag = dose_flag, n.chains = n.chains,
                                       inits.list.set = inits.list.set)
          if(n.chains > 1){
            ## use the 'potential scale reduction factor' to check convergence mixing, reference
            ## https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
            diag.converge <- c()
            for(k in 1:length(retrieve_param)){
              diag.converge <- c(diag.converge,
                                 gelman.diag(as.mcmc.list(post_samples[[k]]),
                                             autoburnin = F)$psrf[ , 1])
            }
            if(any(diag.converge > 1.1)){
              message <- paste(message, ("\n MCMC fail to converge"))}
          }

          ## dose recommendation ##
          if(dose_flag == 1){

            ## if only one dose is assigned in the study, need to modify later

            sim.betas <- as.matrix(rbind(post_samples$beta_other[, , 1],
                                         post_samples$beta_dose[, , 1]))
            mnTTP.dose1 <-
              mean(apply(sim.betas, 2,
                         function(o) { as.numeric(o[1] + o[2] < thrd1)}))
            ## mnTTP.dose1: Pr[(dose = 1) + (cycle = 1) < 0.28]

            if(mnTTP.dose1 < p_tox1) {
              message <- paste(message,
                               ("\n early stop, the recommended dose is 1"))
              break.label = T
              break
            } else {
              doseA <- doseA + 1    ## increase the dose
              dose_flag = 0

              ### Don't change dose level in the next cycle
              rec_dose_act <- rec_dose_act  ##
              cycle_act <- cycle_act + 1    ## update cycle
            }
          } else {

            ## safe dose determination, for early termination##
            if(n.cohort > 2){
              # The termination only works for cohort >= 3 #
              allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                           doses = doses, cycles = cycles,
                                           thrd1 = thrd1, thrd2 = thrd2,
                                           p_tox1 = p_tox1, p_tox2 = p_tox2,
                                           ICD.flag = ICD.flag)
              if (length(allow.doses) == 0) {
                message <- paste(message, "\n early stop, no allowable dose level")
                break.label = T
                break
              }
            }

            Max_tested_doseA = max(rec_doseA) ## calculate the maximum tested dose level for cycle1

            ## recommend dose for new cohort ##

            doseA <- stg1.dos.rec(post_samples = post_samples, doses = doses,
                                  tox.target = tox.target,
                                  Max_tested_doseA = Max_tested_doseA)

            ################################################
            # update the active patient for the next cycle #
            ################################################

            if(ICD.flag == T){

              ### allow dose modification ###

              ## if more than one dose is assigned in the study
              uniq_ID <- unique(patlist$PatID)

              ## safe dose for active patients for next cycle ##
              dos.rec.i.result <-
                stg1.dos.rec.i(post_samples = post_samples, uniq_ID = uniq_ID,
                               patID_act = patID_act, cycle_act = cycle_act,
                               rec_dose_act = rec_dose_act,
                               Max_tested_doseA = Max_tested_doseA,
                               doses = doses, cycles = cycles, c1 = tox.target,
                               p1 = ICD_thrd, DLT.drop.flag = DLT.drop.flag,
                               y.dlt = outcome["y.dlt", ], testedD = testedD)

              cycle_act <- dos.rec.i.result$cycle_nxt
              patID_act <- dos.rec.i.result$patID_nxt
              rec_dose_act <- dos.rec.i.result$rec_dose_nxt
            } else{

              ### no dose modification ###

              cycle_act <- cycle_act + 1
              if(DLT.drop.flag == T){
                act.index <- which(cycle_act <= 6 & outcome["y.dlt", ] == 0)
              } else {
                act.index <- which(cycle_act <= 6)
              }
              cycle_act <- cycle_act[act.index]
              patID_act <- patID_act[act.index]
              rec_dose_act <- rec_dose_act[act.index]
            }
          }
        } else {
          ### no dose modification ###
          cycle_act <- cycle_act + 1
          if(DLT.drop.flag == T){
            act.index <- which(cycle_act <= 6 & outcome["y.dlt", ] == 0)
          } else {
            act.index <- which(cycle_act <= 6)
          }
          cycle_act <- cycle_act[act.index]
          patID_act <- patID_act[act.index]
          rec_dose_act <- rec_dose_act[act.index]
        }
      } else {
        ### stage 2 ###
        if(n.cohort < Max.cohort | IED.flag == T |
           (ICD.flag == T & IED.flag == F)){

          ## model fitting (Skip model fitting when ICD.flag turn off) ##
          post_samples <- phase1stage2(patlist = patlist,
                                       ctrl_param = ctrl_param,
                                       n.iters = n.iters - burn.in,
                                       burn.in = burn.in,
                                       retrieve_param =
                                         c("beta_dose", "beta_other", "alpha", "gamma"),
                                       n.chains = n.chains,
                                       dose_flag = dose_flag)
          ## safe dose determination
          if(ICD.flag == T & IED.flag == F){
            # turn off the ICD.flag option in this scenario
            allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                         doses = doses, cycles = cycles,
                                         thrd1 = thrd1, thrd2 = thrd2,
                                         p_tox1 = p_tox1, p_tox2 = p_tox2,
                                         ICD.flag = F)
          }else{
            allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                         doses = doses, cycles = cycles,
                                         thrd1 = thrd1, thrd2 = thrd2,
                                         p_tox1 = p_tox1, p_tox2 = p_tox2,
                                         ICD.flag = ICD.flag)
          }

          if (length(allow.doses) == 0) {
            message <- paste(message, "\n early stop, no allowable dose level")
            break.label = T
            break
          }

          Max_tested_doseA = max(rec_doseA) ## calculate the maximum tested dose level for cycle1

          ## dose recommendation for new cohort##
          doseA <- stg2.eff.rec(post_samples = post_samples,
                                allow.doses = allow.doses,
                                Max_tested_doseA = Max_tested_doseA,
                                proxy.thrd = proxy.thrd)
        }


        ################################################
        # update the active patient for the next cycle #??? IED fomulation...
        ################################################

        if(IED.flag == T){

          ### allow dose modification ###

          ## if more than one dose is assigned in the study
          uniq_ID <- unique(patlist$PatID)

          ## safe dose for active patients for next cycle ##
          dos.rec.i.result <-
            stg1.dos.rec.i(post_samples = post_samples, uniq_ID = uniq_ID,
                           patID_act = patID_act, cycle_act = cycle_act,
                           rec_dose_act = rec_dose_act,
                           Max_tested_doseA = Max_tested_doseA,
                           doses = doses, cycles = cycles, c1 = tox.target,
                           p1 = ICD_thrd, DLT.drop.flag = DLT.drop.flag,
                           y.dlt = outcome["y.dlt", ], testedD = testedD)

          cycle_act <- dos.rec.i.result$cycle_nxt
          patID_act <- dos.rec.i.result$patID_nxt
          dos.rec.i.IED <- stg2.eff.rec.i(
            post_samples, cycle_nxt = dos.rec.i.result$cycle_nxt,
            rec_dose_nxt = dos.rec.i.result$rec_dose_nxt,
            proxy.thrd = proxy.thrd)
          rec_dose_act <- dos.rec.i.IED$rec_dose_nxt
        } else {

          cycle_act <- cycle_act + 1
          if(DLT.drop.flag == T){
            act.index <- which(cycle_act <= 6 & outcome["y.dlt", ] == 0)
          } else {
            act.index <- which(cycle_act <= 6)
          }
          cycle_act <- cycle_act[act.index]
          patID_act <- patID_act[act.index]

          if(ICD.flag == T){
            recom.MED.cycle <- stg2.eff.rec.i.all.cycle(
              post_samples = post_samples, allow.doses = allow.doses,
              cycles = cycles, proxy.thrd = proxy.thrd)$recom.MED.cycle
            rec_dose_act <- recom.MED.cycle[cycle_act]
          }else{
            ### no dose modification ###
            rec_dose_act <- rec_dose_act[act.index]
          }
        }
      }

      if(n.cohort <= Max.cohort){n.cohort <- n.cohort + 1}     ## increase the number of cohorts
    }

    #############################
    # model fitting for stage 3 #
    #############################
    if(break.label == F){
      # study dosen't break before stage3
      ### stage 3 ###
      ## model fitting ##

      if(effcy.flag == T){
        retrieve_param = c("beta_dose", "beta_other", "alpha", "gamma")
        post_samples <- phase1stage2(patlist = patlist,
                                     ctrl_param = ctrl_param,
                                     n.iters = n.iters - burn.in,
                                     burn.in = burn.in,
                                     retrieve_param = retrieve_param,
                                     n.chains = ifelse(n.chains == 1, 2, n.chains),
                                     # need at least two MCMC chains to check the convergence of MCMC chain
                                     dose_flag = dose_flag)

        ## Check MCMC convergency ##
        diag.converge <- c()
        for(k in 1:length(retrieve_param)){
          diag.converge <- c(diag.converge,
                             gelman.diag(as.mcmc.list(post_samples[[k]]),
                                         autoburnin = F)$psrf[ , 1])
        }
        if(any(diag.converge > 1.1)){
          message <- paste(message, ("\n MCMC fail to converge"))}

        ## safe dose determination
        allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                     doses = doses, cycles = cycles,
                                     thrd1 = thrd1, thrd2 = thrd2,
                                     p_tox1 = p_tox1, p_tox2 = p_tox2,
                                     ICD.flag = ICD.flag)

        if (length(allow.doses) == 0) {
          message <- paste(message, "\n early stop, no allowable dose level")
          break.label = T
        } else {
          ##################
          # dose recommend #
          ##################
          Max_tested_doseA = max(rec_doseA) ## calculate the maximum tested dose level for cycle1
          doseA <- stg3.eff.rec(post_samples = post_samples,
                                allow.doses = allow.doses,
                                Max_tested_doseA = Max_tested_doseA,
                                proxy.thrd = proxy.thrd)
        }
      } else {
        ## fit model without efficacy data ##
        retrieve_param <- c("beta_dose", "beta_other", "gamma")
        post_samples <- phase1stage1(patlist = patlist,
                                     ctrl_param = ctrl_param,
                                     n.iters = n.iters - burn.in,
                                     burn.in = burn.in,
                                     retrieve_param = retrieve_param,
                                     dose_flag = dose_flag,
                                     n.chains = ifelse(n.chains == 1, 2, n.chains),
                                     inits.list.set = inits.list.set)
        diag.converge <- c()
        for(k in 1:length(retrieve_param)){
          diag.converge <- c(diag.converge,
                             gelman.diag(as.mcmc.list(post_samples[[k]]),
                                         autoburnin = F)$psrf[ , 1])
        }
        if(any(diag.converge > 1.1)){
          message <- paste(message, ("\n MCMC fail to converge"))}

        ## safe dose determination
        allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                     doses = doses, cycles = cycles,
                                     thrd1 = thrd1, thrd2 = thrd2,
                                     p_tox1 = p_tox1, p_tox2 = p_tox2,
                                     ICD.flag = ICD.flag)

        if (length(allow.doses) == 0) {
          message <- paste(message, "\n early stop, no allowable dose level")
          break.label = T
        } else {

          ##################
          # dose recommend #
          ##################

          Max_tested_doseA = max(rec_doseA) ## calculate the maximum tested dose level for cycle1

          ## recommend dose for new cohort ##

          doseA <- stg3.dos.rec(post_samples = post_samples, doses = doses,
                                tox.target = tox.target,
                                Max_tested_doseA = Max_tested_doseA)
        }
      }

      pp.nTTPM <- pp.nTTP(post_samples, doses, cycles, tox.target)

      ## save the simulation result, dose alocation result and recommended dose level
      patlist$effz <- NULL
      list_simul[[i]] <- list(patlist = patlist, dose_aloca = dose_aloca,
                              doseA = ifelse(break.label, "early break", doseA),
                              n.cohort = Max.cohort, pp.nTTPM = pp.nTTPM,
                              message = message)
      # don't recommmend dose if no allowable dose from stage3 model fitting
    } else {
      ## If the study break before entering stage 3, just save the simulation data
      list_simul[[i]] <- list(patlist = patlist, dose_aloca = dose_aloca,
                              doseA = "early break",
                              n.cohort = ifelse(n.cohort < Max.cohort,
                                                n.cohort, Max.cohort),
                              message = message)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  sim.time <- proc.time() - t
  res <- list(senerio_sum = senerio_sum, eff_sum = eff_sum,
              list_simul = list_simul, doses = doses, cycles = cycles,
              chSize = chSize, sim.time = sim.time, effcy.flag = effcy.flag,
              proxy.thrd = proxy.thrd, DLT.drop.flag = DLT.drop.flag)
  attr(res,'class') <- 'SimPRMD'
  return(res)
}



