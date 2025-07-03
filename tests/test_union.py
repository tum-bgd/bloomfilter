import bgdbloomfilter as bf
import numpy as np
import random
import string

Ncol = 128
N=1024*Ncol


###
table1 = np.random.randint(-100000,100000,size=N).reshape(-1,Ncol)
table2 = np.random.randint(-100000,100000,size=N).reshape(-1,Ncol)
table3 = np.random.randint(-100000,100000,size=N).reshape(-1,Ncol)

## Create a Bloom Filter
f = bf.bloomfilter()
g = bf.bloomfilter()
## compute for certain FPR
fpr=0.01
n = len(table1)
m = int(-n*np.log(fpr) / (np.log(2)*np.log(2)));
k = int(m/n *np.log(2))
print("fpr=%.2f, n=%.2f,m=%.2f, k=%.2f"%(fpr, m,n,k))


f.configure(k,m)
f.insert(table1)
g.configure(k,m)
g.insert(table2)

f.union(g);


print("Regression: All inserted is: %s " % np.all(f.test(table2)))
