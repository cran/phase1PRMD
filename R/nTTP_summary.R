
nTTP_summary <- function(Tox.prob.M, nTTP.all, wm) {

  #' Generate the mean nTTP score and the probability of observing DLT for all
  #' doses and cycles
  #'
  #' \code{nTTP_summary} generates the mean nTTP score and the probability of
  #' observing DLT for all doses and cycles
  #'
  #' @param Tox.prob.M      Toxicity probability matrix with 4 dimension: dose,
  #'   cycle, type, grade. Tox.prob.M can be the output of the build-in matrix
  #'   of function \code{\link[phase1RMD]{GenToxProb}} in package phase1RMD.
  #'   See more details about how to generate toxicity probability matrices in
  #'   the help document of \code{\link[phase1RMD]{GenToxProb}}.
  #' @param nTTP.all The output of \code{\link{nTTP.array}}
  #' @param wm      (numerical matrix) Toxicity weighted matrix, with row be the
  #'   type of the toxicity and column be the toxicity grade
  #'
  #' @return \item{mnTTP.M}{matrix of mean nTTP for all doses and cycles}
  #'   \item{pDLT.M}{matrix of probability of observing DLT for all doses and
  #'   cycles}
  #'
  #' @examples
  #'
  #' data("prob")
  #'
  #' wm <- matrix(c(0, 0.5, 0.75, 1, 1.5,
  #'                0, 0.5, 0.75, 1, 1.5,
  #'                0, 0, 0, 0.5, 1),
  #'              byrow = TRUE, ncol = 5)          # weighted matrix for toxicity matrix
  #'                                               # nrow = No.of type; ncol = No. of grade
  #' toxmax <- 2.5
  #'
  #' nTTP.all <- nTTP.array(wm, toxmax)
  #'
  #' tox.matrix <- prob["MTD4", "flat", , , , ]
  #'
  #' nTTP_summary(tox.matrix, nTTP.all, wm)
  #'
  #' @importFrom arrayhelpers vec2array
  #' @export

  ## precalculation
  dims <- dim(Tox.prob.M)
  n.dose <- dims[1]
  n.cycle <- dims[2]
  n.tox.type <- dims[3]
  n.tox.grade <- dims[4]

  ## preallocation
  mnTTP.M <- matrix(NA, n.dose, n.cycle)
  pDLT.M <- matrix(NA, n.dose, n.cycle)

  nTTProb <- array(NA, c(rep(n.tox.grade, n.tox.type)))
  capacity <- n.tox.grade^n.tox.type

  for(i in 1:n.dose){
    for(j in 1:n.cycle){

      for(k in 1:capacity){
        index <- vec2array(k, dim = c(rep(n.tox.grade, n.tox.type)))
        nTTProb[k] <- prod(Tox.prob.M[i, j, ,][cbind(1:n.tox.type, c(index))])
      }

      mnTTP.M[i, j] <- mnTTP(nTTProb, nTTP.all)
      pDLT.M[i, j] <- pDLT(nTTProb, wm)
    }
  }

  rownames(mnTTP.M) <- paste0("D", 1:n.dose)
  rownames(pDLT.M) <- paste0("D", 1:n.dose)
  colnames(mnTTP.M) <- paste0("C", 1:n.cycle)
  colnames(pDLT.M) <- paste0("C", 1:n.cycle)

  return(list(mnTTP.M = mnTTP.M, pDLT.M = pDLT.M))
}
