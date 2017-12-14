# `sam.fit`: Fit SAM model

## Description


 Fit SAM model


## Usage

```r
sam.fit(data, conf, parameters, newtonsteps = 3, rm.unidentified = FALSE,
  run = TRUE, lower = getLowerBounds(parameters),
  upper = getUpperBounds(parameters), sim.condRE = TRUE, ...)
```


## Arguments

Argument      |Description
------------- |----------------
```data```     |     data for the sam model as returned from the setup.sam.data function
```conf```     |     model configuration which can be set up using the [`defcon`](defcon.html) function and then modified either directly in R or by saving it to a text file using the function [`saveConf`](saveConf.html) , modifying the text file, and then reading the configuration from the textfile using the function [`loadConf`](loadConf.html) . For more details about the configuration see details.
```parameters```     |     initial values which can be set up using the [`defpar`](defpar.html) function and then modified.
```newtonsteps```     |     optional extra true newton steps
```rm.unidentified```     |     option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
```run```     |     if FALSE return AD object without running the optimization
```lower```     |     named list with lower bounds for optimization (only met before extra newton steps)
```upper```     |     named list with upper bounds for optimization (only met before extra newton steps)
```sim.condRE```     |     logical with default `TRUE` . Simulated observations will be conditional on estimated values of F and N, rather than also simulating F and N forward from their initial values.
```...```     |     extra arguments to MakeADFun

## Details


 The model configuration object `conf` is a list of different objects defining different parts of the model. The different elements of the list are:
 list("\n", list(list("$minAge:"), list("A single integer defining the the lowest age class in the assessment.")), "\n", list(list("$maxAge:"), list("A single integer defining the the highest age class in the assessment.")), "\n", list(list("$maxAgePlusGroup:"), list("Is last age group considered a plus group (1 yes, or 0 no).")), "\n", list(list("$keyLogFsta:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of the fishing mortality states (normally only first row is used). '-1' is used for entries where no fishing mortality applies (e.g. age groups in survey fleets, or unobserved age groups). For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding state vector. If the same number is used for two age classes, then the fishing mortality for those age classes are assumed equal (linked to the same state).")), 
    "\n", list(list("$corFlag:"), list("A single integer to specify the correlation structure of log-scale fishing mortality increments (0 independent, 1 compound symmetry, or 2 AR(1)).")), "\n", list(list("$keyLogFpar:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of survey catchability parameters (so only used for survey fleets). '-1' is used for entries where catchability should not be specified (e.g. fleet - age groups combinations where fishing mortality is specified above, or unobserved fleet - age group combinations).  For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal (linked to the same parameter).")), 
    "\n", list(list("$keyQpow:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of density dependent catchability power parameters. This can only be applied to fleets - age combinations where a catchability is defined. '-1' is used for entries where this cannot be applied (e.g. fleet - age groups combinations where fishing mortality is specified above, or unobserved fleet - age group combinations). '-1' is also used to specify that density dependent catchability power parameters is turned off (the most common setup). For entries where density dependent catchability power parameter is to be estimates entries consecutive integers starting at zero must be used. If the same number is used for two age classes, then the density dependent catchability power parameter for those age classes are assumed equal (linked to the same parameter).")), 
    "\n", list(list("$keyVarF:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of variance parameters for the different states in the log-scale fishing mortality random walk process. '-1' should be used for entries where no fishing mortality state is defined in ", list("keyLogFsta"), " above. For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal (linked to the same parameter). ((a curiosity of this setup is that it is possible to set different variance parameter indices for F-states that are coupled in ", 
        list("keyLogFsta"), ". This is ignored and the index corresponding to the lowest F-state number is used)).")), "\n", list(list("$keyVarLogN:"), list("A vector of integers. The length of the vector is equal to the number of age classes. The vector describes the coupling of variance parameters for the log(N)-process. Consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal. A typical setup is to use a unique index for the first age group, because that corresponds to the variance in the (stock-)recruitment, which is often not similar to the variance in the survival process from year to year.")), 
    "\n", list(list("$keyVarObs:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of observation variance parameters. '-1' should be used for entries where no observations are available. For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the observation variance for those age classes are assumed equal (linked to the same parameter).")), 
    "\n", list(list("$obsCorStruct:"), list("A factor specifying the covariance structure used across ages for each fleet. The length of the factor is equal to the number of fleets. The possible options are: (\"ID\" independent, \"AR\" AR(1), or \"US\" for unstructured).")), "\n", list(list("$keyCorObs:"), list("A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes _minus_ _one_. The matrix describes the coupling AR correlations between age classes, and hence is only meaningful for fleets where the \"AR\" observation correlation structure is chosen. '-1' should be used for entries where no observations are available. Notice that the matrix has one column less than the number of age classes, which is because the correlation between age classes is described. Consecutive integers starting at zero must be used. If the same number is used for a given fleet it means that a normal AR(1) structure is used. If different numbers are used for a fleet it means that the correlation parameter changes where the numbers differ. If the \"AR\" structure is specified above, then the corresponding row in this matrix must have valid non-negative entries.")), 
    "\n", list(list("$stockRecruitmentModelCode:"), list("A single integer to specify the stock recruitment connection to use (0 for plain random walk on log recruitment, 1 for Ricker, and 2 for Beverton-Holt).")), "\n", list(list("$noScaledYears:"), list("A single integer specifying the number of years where catch scaling is to be estimated (most often 0, as this is a somewhat exotic option).")), "\n", list(list("$keyScaledYears:"), list("A vector of the years where catch scaling is applied (length should match ", 
        list("noScaledYears"), ") (most often empty, as this is a somewhat exotic option).")), "\n", list(list("$keyParScaledYA:"), list("A matrix of integers specifying the couplings of scale parameters (nrow = ", list("noScaledYears"), ", ncols = no ages) (most often empty, as this is a somewhat exotic option).")), "\n", list(list("$fbarRange:"), list("An integer vector of length 2 specifying lowest and highest age included in Fbar (average fishing mortality summary).")), "\n", list(list("$keyBiomassTreat:"), 
        list("A vector of integers with length equal to the number of fleets. '-1' should be used for entries where the corresponding fleet is not a mass index. A the corresponding fleet is a mass index, then three options are available (0 SSB index, 1 catch index, and 2 FSB index).")), "\n", list(list("$obsLikelihoodFlag:"), list("A factor specifying the type of likelihood to use for each fleet. The length of the factor is equal to the number of fleets. The possible options are: (\"LN\" for log-normal and \"ALN\" Additive logistic normal).")), 
    "\n", list(list("$fixVarToWeight:"), list("A single integer. If weight attribute is supplied for observations this option defines how it is treated (0 as relative weight, 1 as a fixed variance = weight).")), "\n") 


## Value


 an object of class `sam` 


## Examples

```r 
 data(nscodData)
 data(nscodConf)
 data(nscodParameters)
 fit <- sam.fit(nscodData, nscodConf, nscodParameters)
 ``` 

