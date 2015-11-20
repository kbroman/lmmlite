#!/usr/bin/env python

import lmm
import numpy

kinship = numpy.genfromtxt('kinship.csv', delimiter=',')
pheno =   numpy.genfromtxt('pheno.csv', delimiter=',')
covar =   numpy.genfromtxt('covar.csv', delimiter=',')

def print_result(result, index, method):
    result_arr = [result[0], result[1][0,0], result[1][1,0], result[2][0,0], result[3]]
    for j in range(len(result_arr)):
        result_arr[j] = str(result_arr[j])
    print(','.join([str(index), method] + result_arr))

print('index,method,hsq,intercept,sex,sigmasq,loglik') # header
for i in range(pheno.shape[1]):
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=True)
    print_result(result, i+1, 'reml')
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=False)
    print_result(result, i+1, 'ml')
