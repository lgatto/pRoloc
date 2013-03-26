# Calculus of the perturbations for one label and for each individuals
# Called from testStep function
# 
# @title loopInTestStep
# @param A a matrix
# @param B the kernel matrix
# @param ind an integer that corresponds to the number of individuals in the testSet
# @return a column-vector that corresponds to the perturbation for each of 
#         the 'ind' individuals
# @author Samuel Wieczorek
# 
loopInTestStep <- function (A, B, ind) 
.Call("file25fd75026a26", A, B, ind, PACKAGE = "pRoloc")

