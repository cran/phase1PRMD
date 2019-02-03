# phase1PRMD

Implements Bayesian phase I repeated measurement design that accounts for multidimensional toxicity endpoints and 
longitudinal efficacy measure from multiple treatment cycles and allows individualized dose modification. The package 
provides a number of model-based phase I designs, including 1 stage 
models with or without individualized dose modification, 3-stage models with or without
individualized dose modification, etc. Functions are provided to recommend
dosage selection based on the data collected in the available patient cohorts
and to simulate trial characteristics given design parameters.

## Installation
This package can be installed directly from R by running the following code:
```
library(devtools)
install_github("LuZhangstat/phase1PRMD")
```

## Usage
```
library(phase1PRMD)
?RunPRMD
?SimPRMD
```
Detailed manual will relase in the near future
