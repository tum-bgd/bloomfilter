import bloomfilter as bf
import numpy as np
import random
import string
Ncol=8
N=1024*1024


###
table1 = np.random.rand(N,Ncol)
table2 = np.random.rand(N,Ncol)

## Create a Bloom Filter
f = bf.bloomfilter()

## compute for certain FPR
fpr=0.01
n = len(table1)
m = int(-n*np.log(fpr) / (np.log(2)*np.log(2)));
k = int(m/n *np.log(2))
print("fpr=%.2f, n=%.2f,m=%.2f, k=%.2f"%(fpr, m,n,k))


f.configure(k,m)
f.insert(table1)
print("Regression: All inserted is: %s " % np.all(f.test(table1)))
print("Random Table Test should approach FP rate: %f " % np.mean(f.test(table2)))
