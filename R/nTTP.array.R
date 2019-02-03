nTTP.array <- function(wm, toxmax) {

  #' Generate the nTTP dictionary
  #'
  #' \code{nTTP.array} generates the nTTP dictionary for all combination of
  #' toxicity type and grade with a given toxicity weighted matrix. Used in
  #' function \code{\link{nTTP_summary}} for checking the toxicity scenario.
  #'
  #' @param wm      (numeric matrix, m by n) Toxicity weighted matrix, with row
  #'   be the type of the toxicity and column be the toxicity grade
  #' @param toxmax  (scalar)         Normalized constant for nTTP
  #'
  #' @return An m dimensional array with dimension \eqn{(n, n, \ldots, n)}. The
  #'   \eqn{(d1, d2 ,\ldots,dm), {di, i = 1\ldots, m} \in (1, \ldots, n)}th
  #'   element is the nTTP when the grade of \eqn{i}th type of toxicity has
  #'   \eqn{di}th toxicity grade.
  #'
  #' @examples
  #' wm = matrix(c(0, 0.5, 0.75, 1, 1.5,
  #'               0, 0.5, 0.75, 1, 1.5,
  #'               0, 0, 0, 0.5, 1),
  #'             byrow = TRUE, ncol = 5)          # weighted matrix for toxicity matrix
  #'                                              # nrow = No.of type; ncol = No. of grade
  #' toxmax = 2.5
  #'
  #' nTTP.array(wm, toxmax)
  #'
  #' @importFrom arrayhelpers vec2array
  #' @export


  n.tox.grade <- ncol(wm)  # retrieve the No. of tox grade
  n.tox.type  <- nrow(wm)  # retrieve the NO. of tox type

  nTTP <- array(NA, c(rep(n.tox.grade, n.tox.type)))
  capacity <- n.tox.grade^n.tox.type
  for(tt in 1 : capacity) {
    index <- vec2array(tt, dim = c(rep(n.tox.grade, n.tox.type)))
    nTTP[tt] <- sqrt(sum(wm[t(rbind(1 : n.tox.type, index))]^2)) / toxmax
  }
  return(nTTP)
}
