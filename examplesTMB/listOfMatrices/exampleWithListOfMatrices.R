library(TMB)
library(Matrix)

compile("exampleWithListOfMatrices.cpp")
dyn.load(dynlib("exampleWithListOfMatrices"))

M1 = matrix(1:4,2,2)
M2 = matrix(1:9,3,3)


data <- list()
data$x <- list(a=M1, b=M2)

par = list(
  test = 0#Not used, but need at leas one parameter it seems.
)
obj <- MakeADFun(data,parameters = par,DLL = "exampleWithListOfMatrices")

