## summary and plot functions for RunPRMD ##

summary.RunPRMD <- function(object, ...){

  #' Summary a RunPRMD object
  #'
  #' Summary a RunPRMD object. Print the information of recommended dosage
  #' selection along with the mean nTTP and the number of DLT for all doses and
  #' cycles. Will print the mean efficacy for all doses and cycles when
  #' implementing \code{\link{RunPRMD}} with option \code{effcy.flag = TRUE}. The
  #' collected data is displayed in a human-readable table whose cell contain 3
  #' values including observed nTTP, DLT, and dose assignment. The higher the
  #' dose, the warmer the cell background color is. The black color of the
  #' records indicates DLT equals 1.
  #'
  #' @param object RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @return \item{object}{The output of function \code{RunPRMD}}
  #'   \item{mnttp.M}{The mean nTTP for all doses and cycles}
  #'   \item{dlt.count.M}{The number of DLT for all doses and cycles}
  #'   \item{eff.M}{The mean efficacy for all doses and cycles. Return
  #'   \code{NULL} when \code{object$effcy.flag == TRUE}}
  #'
  #' @export
  #'
  #' @examples
  #' ## Check ?RunPRMD for example
  #'

  mnttp.M <- sapply(object$cycles, function(c){
    sapply(object$doses, function(d){
      mean(object$patlist$nTTP[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
    })
  })
  colnames(mnttp.M) <- paste0("C", object$cycles)
  rownames(mnttp.M) <- paste0("D", object$doses)

  dlt.count.M <- sapply(object$cycles, function(c){
    sapply(object$doses, function(d){
      sum(object$patlist$dlt[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
    })
  })
  colnames(dlt.count.M) <- paste0("C", object$cycles)
  rownames(dlt.count.M) <- paste0("D", object$doses)

  if(object$effcy.flag == T){
    eff.M <- sapply(object$cycles, function(c){
      sapply(object$doses, function(d){
        sum(object$patlist$efficacy[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
      })
    })
  }else{eff.M = NULL}
  ans <- list(object = object, mnttp.M = mnttp.M,
              dlt.count.M = dlt.count.M, eff.M = eff.M)
  class(ans) <- "summary.RunPRMD"
  ans
}

print.summary.RunPRMD <- function(x, ...){

  #' Displays a useful description of a summary.RunPRMD object
  #'
  #' Displays a useful description of a summary.RunPRMD object. Call by
  #' \code{link{summary.RunPRMD}}. Check \code{link{summary.RunPRMD}} for the
  #' details of the print information.
  #'
  #' @param x summary.RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @importFrom RColorBrewer brewer.pal
  #' @importFrom dplyr mutate rowwise mutate_at starts_with select funs vars
  #' @import kableExtra
  #' @import knitr
  #' @importFrom utils capture.output
  #' @export

  . <- NULL; PatID <- NULL # for suppressing the note in checking

  if (!inherits(x, "summary.RunPRMD"))
    stop(gettextf("'x' must inherit from class %s",
                  dQuote("summary.RunPRMD")),
         domain = NA)

  n.cycle <- max(x$object$cycles); n.dose <- max(x$object$doses)
  color.pal <- brewer.pal(n = n.dose, name = "RdYlGn")[n.dose:1]
  l <- length(x$object$patlist$PatID)
  uniq_ID <- unique(x$object$patlist$PatID)

  n.patient <- length(uniq_ID)
  report <- matrix(NA, nrow = n.patient, ncol = n.cycle)
  colnames(report) <- c(paste0("cycle", 1:n.cycle))
  rownames(report) <- uniq_ID

  for(i in 1:l){
    report[which(uniq_ID == x$object$patlist$PatID[i]),
           x$object$patlist$cycle[i]] <-
      paste(format(round(x$object$patlist$nTTP[i],3),3), x$object$patlist$dlt[i],
            x$object$patlist$dose[i], sep = ",")
  }

  report.data <- data.frame(report)

  cat("\n The mean nTTP for all doses and cycles: \n")
  print(format(round(x$mnttp.M, 3), 3))
  cat("\n The number of DLT for all doses and cycles: \n")
  print(format(round(x$dlt.count.M, 3), 3))
  if(x$object$effcy.flag == T){
    cat("\n The mean efficacy for all doses and cycles: \n")
    print(format(round(x$eff.M, 3), 3))
  }
  cat("Recommended dose", x$object$doseA, "\n")
  cat("\nFor patients: \n", x$object$pat_rec$patID,
      "\non cycle: \n", x$object$pat_rec$cycle,
      "\nWe suggest dose levels: \n", x$object$pat_rec$dose, "\n")

  invisible(x)
  invisible(capture.output(print(
  report.data %>%
  mutate(PatID = row.names(.))%>%
  rowwise()%>%
  mutate_at(vars(starts_with("cycle")), funs(
    ifelse(!is.na(.),
           ifelse(as.numeric(unlist(strsplit(as.character(.), ","))[2]) == 0,
                  cell_spec(., color = "white", bold = T,
                            background = color.pal[as.numeric(strsplit(as.character(.), ",")[[1]])[3]]),
                  cell_spec(., color = "black", bold = T,
                            background = color.pal[as.numeric(strsplit(as.character(.), ",")[[1]])[3]])),
           cell_spec(., color = "white", background = "white"))
  ))%>%
  select(PatID, starts_with("cycle"))%>%
  kable(escape = F, format = "html") %>%
  kable_styling() %>%
  footnote((general = "the entries in each cell are formatted in: nTTP, dlt, dose"))
  )))
}

plot.RunPRMD <- function(x, ..., select_cycle = x$cycles){

  #' Plot nTTP and efficacy boxplots of a RunPRMD object
  #'
  #' Plot nTTP boxplots of a RunPRMD object. Plot efficacy boxplots when
  #' implementing \code{RunPRMD} with option \code{effcy.flag == TRUE}.
  #'
  #' @param x RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #' @param select_cycle A vector indication the cycle in the boxplot. Default
  #'   is \code{cycle} of \code{x}.
  #'
  #' @examples
  #' ## Check ?RunPRMD for example
  #'
  #' @import ggplot2
  #'
  #' @export

  dose <- NA; cycle <- NA; nTTP <- NA; effcy <- NA
  thislist <- list(...)
  Max.nTTP <- max(x$patlist$nTTP)
  select.index <- sapply(x$patlist$cycle, function(a){any(a == select_cycle)})
  nttp.dt <- data.frame(
    nTTP = as.vector(x$patlist$nTTP[select.index]),
    cycle = factor(x$patlist$cycle[select.index], levels = select_cycle),
    dose = factor(x$patlist$dose[select.index], levels = x$doses))
  pnTTP <- ggplot(data = nttp.dt, aes(x = dose, y = nTTP,
                                      fill = cycle)) +
    stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() +
    facet_wrap(~cycle) + xlab("Doses") + ylab("nTTP")
    scale_x_discrete(breaks = x$doses, labels = x$doses)
  print(pnTTP)
  if(x$effcy.flag == T){
    readline(prompt="Press [enter] to check the efficacy box-plot")
    effcy.dt <- data.frame(
      effcy = as.vector(x$patlist$efficacy[select.index]),
      cycle = factor(x$patlist$cycle[select.index], levels = select_cycle),
      dose = factor(x$patlist$dose[select.index], levels = x$doses))
    peff <- ggplot(data = effcy.dt, aes(x = dose, y = effcy,
                                        fill = cycle)) +
      stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() +
      facet_wrap(~cycle) + xlab("Doses") + ylab("Efficacy")
      scale_x_discrete(breaks = x$doses, labels = x$doses)
    print(peff)
  }
}


patlist.display <- function(patlist, n.dose, n.cycle){

  #' Display patient records
  #'
  #' Display patient records in a human-readable table. Each cell contains 3
  #' values including observed nTTP, DLT, and dose assignment. The higher the
  #' dose, the warmer the cell background color is. The black color of the
  #' records indicates DLT equals 1.
  #'
  #' @param patlist A list of the patient treatment records, which must contains
  #'   the following variables: \describe{ \item{PatID}{denotes the patient ID
  #'   where the elements are specified by cohort and subject number. For
  #'   example, "cohort2subject3" denotes the third subject in the second
  #'   cohort} \item{dose}{records the dose level assigned for each patient
  #'   through the whole treatment} \item{cycle}{shows the treatment cycle
  #'   information of each record} \item{nTTP}{records the corresponding nTTP
  #'   score.} \item{dlt}{indicates whether a DLT event is observed or not?}}
  #' @param n.dose  The number of dose in the study
  #' @param n.cycle The number of cycle in the study
  #'
  #' @importFrom RColorBrewer brewer.pal
  #' @importFrom dplyr mutate rowwise mutate_at starts_with select funs vars
  #' @import kableExtra
  #' @import knitr
  #' @importFrom utils capture.output
  #' @export

  . <- NULL; PatID <- NULL # for suppressing the note

  color.pal <- brewer.pal(n = n.dose, name = "RdYlGn")[n.dose:1]
  l <- length(patlist$PatID)
  uniq_ID <- unique(patlist$PatID)

  n.patient <- length(uniq_ID)
  report <- matrix(NA, nrow = n.patient, ncol = n.cycle)
  colnames(report) <- c(paste0("cycle", 1:n.cycle))
  rownames(report) <- uniq_ID

  for(i in 1:l){
    report[which(uniq_ID == patlist$PatID[i]),
           patlist$cycle[i]] <-
      paste(format(round(patlist$nTTP[i],3),3), patlist$dlt[i],
            patlist$dose[i], sep = ",")
  }

  report.data <- data.frame(report)

  invisible(capture.output(print(
    report.data %>%
      mutate(PatID = row.names(.))%>%
      rowwise()%>%
      mutate_at(vars(starts_with("cycle")), funs(
        ifelse(!is.na(.),
               ifelse(as.numeric(unlist(strsplit(as.character(.), ","))[2]) == 0,
                      cell_spec(., color = "white", bold = T,
                                background = color.pal[
                                  as.numeric(strsplit(as.character(.), ",")[[1]])[3]]),
                      cell_spec(., color = "black", bold = T,
                                background = color.pal[
                                  as.numeric(strsplit(as.character(.), ",")[[1]])[3]])),
               cell_spec(., color = "white", background = "white"))
      ))%>%
      select(PatID, starts_with("cycle"))%>%
      kable(escape = F, format = "html") %>%
      kable_styling() %>%
      footnote((general = "the entries in each cell are formatted in: nTTP, dlt, dose"))
    )))
}

