#rm(list = ls())
GenToxProb <- function(toxtype = c("Neurological", "Renal", "Hematological"), 
                       intercept.alpha = c(2, 3, 4.2, 5.7), 
                       coef.beta = c(-0.2, -0.4, -0.7), 
                       cycle.gamma = 0, 
                       MaxCycle = 6,
                       Numdose = 6) {
        if(length(toxtype) != length(coef.beta))
                stop("The number of toxicity types should match that of betas!")
        
        if(length(intercept.alpha) != 4){
                stop("Exactly four intercepts alpha are needed for grade 0--4 in simulation!")
        } 
        
        if(min(diff(intercept.alpha, lag = 1, differences = 1)) < 0){
                stop("Intercepts alpha for simulation must be in a monotonic increasing order!")
        }
        
        tox.matrix <- array(NA, dim = c(Numdose, MaxCycle, length(toxtype), 5))
        
        for(i in 1:MaxCycle){   
                for(j in 1:Numdose){ 
                        cump <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cumulative probability
                        celp <- matrix(NA, nrow = length(toxtype), ncol = 5)  # cell probability
                        
                        for(k in 1:length(toxtype)){    
                                # proportional odds model
                                logitcp <- intercept.alpha + coef.beta[k] * j +  cycle.gamma * (i - 1)
                                # cumulative probabilities
                                cump[k, 1:4] <- exp(logitcp) / (1 + exp(logitcp)) 
                                cump[k, 5] <- 1
                                # cell probabilities
                                celp[k, ] <- c(cump[k,1], diff(cump[k, ], lag = 1, differences = 1))
                        }
                        
                        tox.matrix[j, i, ,] <- celp
                }
        }
        
        return(tox.matrix)
}

