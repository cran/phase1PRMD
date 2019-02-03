## summary and plot functions for SimPRMD ##

summary.SimPRMD <- function(object, ...){

  #' Summary a SimPRMD object
  #'
  #' @param object SimPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @return
  #' \item{senerio_sum}{Output \code{senerio_sum} of function \code{SimPRMD}}
  #' \item{eff_sum}{Output \code{eff_sum} of function \code{SimPRMD}}
  #' \item{n.trial}{The number of trials in the simulation}
  #' \item{alloc.perc}{The dose allocation percentage for cycle1}
  #' \item{n.stop}{The number of early stop cases}
  #' \item{m.n.pat}{The average number of patient in each trial}
  #' \item{m.dlt.rt}{dlt rate}
  #' \item{c1_dlt.rt}{Cycle 1 dlt rate}
  #' \item{cs_dlt.rt}{Subsequent cycle (cycle > 1) dlt rate}
  #' \item{alloc.perc}{Dose allocation of cycle 1}
  #' \item{sbsq.alloc}{Dose allocation of subsequent cycles (cycle > 1)}
  #' \item{rec.prec}{The percentage of Recommended doses for cycle 1}
  #' \item{effcy.flag}{Argument \code{effcy.flag} of function \code{SimPRMD}}
  #' \item{DLT.drop.flag}{Argument \code{DLT.drop.flag} of function
  #' \code{SimPRMD}}
  #'
  #'
  #' @examples
  #' ## Check ?SimPRMD for example
  #'
  #' @importFrom stats cor median qt quantile sd var cycle
  #' @export

  thislist <- list(...)
  chSize <- object$chSize
  #-------------allocation percentage-------------------#
  #first cycle
  alloc <- sapply(object$list_simul, function(a){
    a$dose_aloca
  })
  alloc.perc <- apply(alloc, 1, function(a){
    sum(a)/sum(alloc)
  })
  #subsequent cycle (cycle > 1)
  dose.sbsq <- unlist(sapply(object$list_simul, function(a){
    a$patlist$dose[which(a$patlist$cycle!= 1)]}))
  sbsq.alloc <- table(factor(dose.sbsq[which(dose.sbsq != "early break")],
                             levels = 1:dim(object$senerio_sum$mnTTP.M)[1])) /
    length(dose.sbsq)

  #------- dose recommendation for cycle 1 & early stop rate --------#
  dose_A_V <- unlist(sapply(object$list_simul, function(a){a$doseA}))
  n.stop <- sum(dose_A_V == "early break")
  n.coh <- unlist(sapply(object$list_simul, function(a){a$n.cohort}))
  m.n.pat <- mean(n.coh) * chSize
  early.ind <- which(dose_A_V == "early break")
  early.stp.sm <- summary(n.coh[early.ind]) # early stop summary

  rec.perc <- table(factor(dose_A_V[which(dose_A_V != "early break")],
                           levels = object$doses)) / (length(dose_A_V) - n.stop)


  #------- mean percentage of DLT rate --------#
  dlt.rt.dt <- unlist(sapply(object$list_simul, function(a){sum(a$patlist$dlt) /
      (a$n.cohort * chSize)}))
  m.dlt.rt <- mean(dlt.rt.dt)
  # cycle 1 and subsequent cycle (cycle > 1) dlt rate
  c1_dlt.rt <- dlt.rt.c1(object$list_simul, chSize)
  cs_dlt.rt <- dlt.rt.subseq(object$list_simul, chSize)

  # mean dose sum per person
  dose.sum.dt <- unlist(sapply(object$list_simul, function(a){sum(a$patlist$dose) /
      (a$n.cohort * chSize)}))
  m.dose.sum <- mean(dose.sum.dt)

  ## distance from the perfect trajectory
  # MTD trajectory
  idel.tj <- apply(object$senerio_sum$mnTTP.M, 2, function(a){
    which.min( ifelse(a <= 0.28, 0.28 - a, 1))
  })
  if(object$effcy.flag == T){
    # when efficacy is considered, we calculate lowest efficacious doses lower than MTD
    for(c in object$cycles){
      idel.tj[c] <- min(which(object$eff_sum$eff.M[1:idel.tj[c], c] >=
                                (max(object$eff_sum$eff.M[1:idel.tj[c], c]) -
                                   object$proxy.thrd)))
    }
  }
  dis.to.idl <- sapply(object$list_simul, function(a){
    sum(abs(a$patlist$dose - idel.tj[a$patlist$cycle])) /
      (a$n.cohort * chSize)
  })
  m.dis.to.idl <- mean(dis.to.idl)

  ## add the summary of median survival duration, mean(1st quartile, 3rd quartile)
  surv_dur <- sapply(object$list_simul, function(a){
    median(rev(a$patlist$cycle)[!duplicated(rev(a$patlist$PatID))])
  })
  sur_sum <- c(mean(surv_dur), quantile(surv_dur, c(0.25, 0.75)))

  ans <- list(senerio_sum = object$senerio_sum, eff_sum = object$eff_sum,
              n.trial = length(object$list_simul),
              n.stop = n.stop, m.n.pat = m.n.pat, m.dlt.rt = m.dlt.rt,
              c1_dlt.rt = c1_dlt.rt, cs_dlt.rt = cs_dlt.rt,
              m.dose.sum = m.dose.sum, m.dis.to.idl = m.dis.to.idl,
              sur_sum = sur_sum, alloc.perc = alloc.perc,
              sbsq.alloc = sbsq.alloc, rec.perc = rec.perc,
              effcy.flag = object$effcy.flag,
              DLT.drop.flag = object$DLT.drop.flag)
  class(ans) <- "summary.SimPRMD"
  ans
}

print.summary.SimPRMD <- function(x,
                                     ...){

  #' Displays a useful description of a summary.SimPRMD object
  #'
  #' Displays a useful description of a summary.SimPRMD object
  #'
  #' @param x summary.SimPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @export

  if (!inherits(x, "summary.SimPRMD"))
    stop(gettextf("'x' must inherit from class %s", dQuote("summary.SimPRMD")),
         domain = NA)

  cat("\n The therotical mean nTTP matrix of all dose levels and cycles are:\n")
  print(x$senerio_sum$mnTTP.M)
  cat("\n The therotical probability of observing DLT of all dose levels and cycles are: \n")
  print(x$senerio_sum$pDLT.M)
  if(x$effcy.flag == TRUE){
    cat("\n The mean efficacy for all dose levels and cycles are: \n")
    print(x$eff_sum$eff.M)
  }
  cat("\nIn total", x$n.trial, "simulations, with", x$n.stop,
      "early stop cases ")
  cat("\nOn average,", x$m.n.pat, "patients are enrolled in each trial")
  if(x$DLT.drop.flag == T){
    cat("\nOn average, the DLT drop off rate is ", round(x$m.dlt.rt*100, 2),
        "%")
    cat("\nDLT drop off rate on Cycle 1 is ", round(x$c1_dlt.rt*100, 2),
        " DLT drop off rate on cycle > 1 is ", round(x$cs_dlt.rt*100, 2))
  }
  cat("\n On average, the culmulative dose per patient is", x$m.dose.sum,
      " doses")
  cat("\n On average, the distance from the perfect trajectory is ",
      x$m.dis.to.idl)
  cat("\n\n The summary of median treatment duration (in cycle): \n")
  print(format(round(x$sur_sum, 3), 3))
  cat("\n The dose allocation for cycle 1 of this simulation is: \n")
  print(paste0(format(round(x$alloc.perc*100, 3), 3), "%"))
  cat("\n The dose allocation for cycle > 1 of this simulation is: \n")
  print(paste0(format(round(x$sbsq.alloc*100, 3), 3), "%"))
  cat("\n The dose recommendation for cycle 1 of this simulation is: \n")
  print(paste0(format(round(x$rec.perc*100, 3), 3), "%"))
  invisible(x)
}

plot.SimPRMD <- function(x, ..., title.add = TRUE){

  #' Plots of a SimPRMD object
  #'
  #' Plot the predictive probability of nTTP < target toxicity for all cycles
  #' and doses , the mean nTTP vs cycle1 and cycle > 2 for all doses of a
  #' SimPRMD object. Plot median treatment duration boxplot along with the DLT
  #' drop off rate when implementing \code{\link{SimPRMD}} with option
  #' \code{DLT.drop.flag = TRUE}.
  #'
  #' @param x SimPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #' @param title.add controls whether there is a title on plots or not.
  #'
  #' @examples
  #' ## Check ?SimPRMD for example
  #'
  #' @import ggplot2
  #' @importFrom gridExtra grid.arrange
  #' @importFrom plyr ddply summarise
  #'
  #' @export

  ## calculate the posterior probability of nttp < target

  cycle <- NA; dose <- NA; ppnttp.d <- NA; q.2.5 <- NA; q.97.5 <- NA
  perc <- NA; surv_dur <- NA; nttp.c <- NA

  dose_A_V <- unlist(sapply(x$list_simul, function(a){a$doseA}))
  n.doses <- length(x$doses)
  n.cycles <- length(x$cycles)
  ppnttp.d <- sapply(which(dose_A_V != "early break"), function(a){
      x$list_simul[[a]]$pp.nTTPM
  })

  ppnttp.dt <- data.frame(
    ppnttp.d = as.vector(ppnttp.d),
    cycle = factor(rep(x$cycles, length(ppnttp.d)/n.cycles),
                   levels = x$cycles),
    dose = factor(rep(rep(x$doses, each = n.cycles),
                      length(ppnttp.d) / (n.doses*n.cycles)),
                  levels = x$doses))
  ppnttp.sm <- summarySE(ppnttp.dt, measurevar="ppnttp.d",
                         groupvars=c("cycle", "dose"))
  pd <- position_dodge(0.6)
  pp <- ggplot(ppnttp.sm, aes(x = cycle, y = ppnttp.d,
                              group = dose,
                              colour = dose)) +
    geom_errorbar(aes(ymin = q.2.5, ymax = q.97.5),
                  width = .4, position = pd) +
    geom_line(position = pd, aes(linetype = dose)) +
    geom_point(position = pd, aes(shape = dose), size = 3,
               fill = "white") +
    scale_y_continuous(limits = c(0.0, 1.0))

  if(title.add == T){
    print(pp +
            ggtitle("predictive probability of nTTP < target toxicity"))
  }else{
    print(pp)
  }
  readline(prompt="Press [enter] to check the next plot")

  ## mean nTTP vs cycle1 and cycle > 2 ##
  m.nttp.c.ls <- list()
  for(cyc in 1:2){
    if(cyc == 1){
      m.nttp.c.ls[[cyc]] <- unlist(sapply(x$list_simul, function(a){
      nttp.dt <- data.frame(nttp.c = a$patlist$nTTP[
        which(a$patlist$cycle == cyc)],
        do.c = a$patlist$dose[which(a$patlist$cycle == cyc)])
      dt <- ddply(nttp.dt, ~do.c, summarise, mean = mean(nttp.c))
      return(dt)
      }))
    } else {
      m.nttp.c.ls[[cyc]] <- unlist(sapply(x$list_simul, function(a){
      nttp.dt <- data.frame(nttp.c = a$patlist$nTTP[
        which(a$patlist$cycle >= cyc)],
        do.c = a$patlist$dose[which(a$patlist$cycle >= cyc)])
      dt <- ddply(nttp.dt, ~do.c, summarise, mean = mean(nttp.c))
      return(dt)
      }))
    }
  }
  m.nttp.dt <- data.frame(m.nttp = unlist(m.nttp.c.ls)[
    which(unlist(m.nttp.c.ls)< 1)],
    dose = factor(unlist(m.nttp.c.ls)[which(unlist(m.nttp.c.ls) >= 1)],
                        levels = x$doses),
    group = c(rep("c 1", length(which(m.nttp.c.ls[[1]] < 1))),
              rep("c >1", length(which(m.nttp.c.ls[[2]] < 1)))))

  tgc <- summarySE(m.nttp.dt, measurevar="m.nttp", groupvars=c("group","dose"))
  pd <- position_dodge(0.4)
  p1 <- ggplot(tgc, aes(x = tgc$dose, y = tgc$m.nttp, group = tgc$group,
                        colour = tgc$group)) +
    geom_errorbar(aes(ymin=tgc$q.2.5, ymax=tgc$q.97.5), width=.4,
                  position=pd) +
    geom_line(position=pd, aes(linetype = tgc$group)) +
    geom_point(position=pd, aes(shape = tgc$group), size = 3, fill = "white") +
    scale_y_continuous(limits = c(0.0, 0.85))

  if(title.add == T){
    print(p1 +
            ggtitle("mean nTTP of all dose levels for cycle 1 vs cycle > 1"))
  }else{
    print(p1)
  }

  if(x$DLT.drop.flag == T){
    # draw the survival plot when DLT.drop.flag is TRUE
    readline(prompt="Press [enter] to check the next plot")
    ## median treatment duration boxplot ##
    surv_dur <- sapply(x$list_simul, function(a){
      median(rev(a$patlist$cycle)[!duplicated(rev(a$patlist$PatID))])
    })
    cycle_table <- table(unlist(sapply(x$list_simul, function(a){
      a$patlist$cycle
    })))
    cycle6_sur_count <- length(unlist(sapply(x$list_simul, function(a){
      a$patlist$dlt[which(a$patlist$dlt == 0 & a$patlist$cycle == 6)]
    })))
    cycle_prec <- c(cycle_table[-1], cycle6_sur_count)/cycle_table[1]

    cycle_data <- data.frame(perc = cycle_prec, cycle = x$cycles)
    surv_data <- data.frame(surv_dur = surv_dur)

    cycle_p_g <- ggplot(cycle_data,
                        aes(x = cycle, y = perc)) +
      geom_bar(stat = "identity", width = 0.3, color = "cadetblue4",
               fill = "azure2", size = 0.5) +
      geom_line(color = "cadetblue4", size = 1) +
      xlab("cycle") + ylab("percentage of safe patients") +
      geom_text(aes(label = paste0(sprintf("%0.2f",
                                           round(100*perc, 2)),
                                   "%")),
                hjust = 0.5, vjust = -1) +
      theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
      scale_x_continuous(minor_breaks = 1, limits = c(0.5, n.cycles + 0.5),
                         breaks = seq(1, n.cycles, 1),
                         labels = seq(1, n.cycles, 1))+
      scale_y_continuous(limits = c(0, 1), breaks =  seq(0, 1, 0.25),
                         labels = paste0(seq(0, 100, 25), "%"))

    surv_dur_box <- ggplot(surv_data, aes(x = "100%", y = surv_dur)) +
      geom_boxplot( color = "cadetblue4", fill = "azure2", size = 1) +
      coord_flip() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin = margin(0, 0, 0, 0, "cm")) +
      xlab("Median") +
      scale_y_continuous(minor_breaks = 1, limits = c(0.5, n.cycles + 0.5),
                         breaks = seq(1, n.cycles, 1),
                         labels = seq(1, n.cycles, 1))
    if(title.add == T){
      grid.arrange(surv_dur_box, cycle_p_g, heights = c(0.9, 3), ncol = 1,
                   nrow = 2, top = "DLT drop off rate and median treatment duration")
    }else{
      grid.arrange(surv_dur_box, cycle_p_g, heights = c(0.9, 3), ncol = 1,
                   nrow = 2)
    }
  }
}


