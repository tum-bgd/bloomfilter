import bloomfilter as bf
import numpy as np
import random
import string
N=1024*1024*8


### Create two large random string tables
def random_string(N):
    return ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(N))

# a random string over a limited alphabet split in one arbitrary character
# provides a nice variable length table of strings
stringtable1 = np.array(random_string(N).split("A")).reshape(-1,1)
stringtable2 = np.array(random_string(N).split("X")).reshape(-1,1)

print(stringtable1)
## Create a Bloom Filter
f = bf.bloomfilter()

## compute for certain FPR
fpr=0.01
n = len(stringtable1)
m = int(-n*np.log(fpr) / (np.log(2)*np.log(2)));
k = int(m/n *np.log(2))
print("fpr=%.2f, n=%.2f,m=%.2f, k=%.2f"%(fpr, m,n,k))



f.configure(k,m)
f.insert(stringtable1)
print("Regression: All inserted is: %s " % np.all(f.test(stringtable1)))

## Get the data
print("Random Table Test should approach FP rate: %f " % np.mean(f.test(stringtable2)))
