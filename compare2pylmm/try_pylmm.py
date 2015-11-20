#!/usr/bin/env python
#
# script to load Recla example data (from R/lmmlite) and do LMM
# calculates with pylmm

import lmm
import numpy

# load the datasets
kinship = numpy.genfromtxt('kinship.csv', delimiter=',')
pheno =   numpy.genfromtxt('pheno.csv', delimiter=',')
covar =   numpy.genfromtxt('covar.csv', delimiter=',')

# function to print results
def print_result(result, index, method):
    result_arr = [result[0], result[1][0,0], result[1][1,0], result[2][0,0], result[3]]
    for j in range(len(result_arr)):
        result_arr[j] = str(result_arr[j])
    print(','.join([str(index), method] + result_arr))

print('index,method,hsq,intercept,sex,sigmasq,loglik') # print header row

# loop over the columns in pheno
#    fit the LMM by REML or ML and print the results to STDIN
for i in range(pheno.shape[1]):
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=True)
    print_result(result, i+1, 'reml')
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=False)
    print_result(result, i+1, 'ml')
