import bloomfilter as bf
import numpy as np

f = bf.bloomfilter();
f.configure(3,1024)


#x = np.array([[1,2],[3,4]])
#print(x)
#f.insert(x)
#print("Strided")
#f.insert(x[:,::2])
#

x = np.array(["Gerne","will","ich", "mich", "bequemen"]).reshape(-1,1)
f.insert(x)

print("INSERT COMPLETE")
for i in range(len(x)):
    print("DBG",f.dbg_getrow(x,i))
    print("DBH",f.dbg_gethex(x,i))

filterdata = f.to_bytes();
print(filterdata)
    
print(f.test(x))
x = np.array(["Kreuz","und","Becher", "anzunehmen"]).reshape(-1,1)
print(f.test(x))
    
g = bf.bloomfilter();
g.from_bytes(filterdata,3,1024);
print("="*63)
print("Testing load filter")
print("="*63)

x = np.array(["Gerne","will","ich", "mich", "bequemen"]).reshape(-1,1)
print(f.test(x))
x = np.array(["Kreuz","und","Becher", "anzunehmen"]).reshape(-1,1)
print(f.test(x))


### Test Pickle support
import pickle
pickle.dump(f, open("test.p","wb"))
g = pickle.load(open("test.p","rb"))


x = np.array(["Gerne","will","ich", "mich", "bequemen"]).reshape(-1,1)
print(g.test(x))
x = np.array(["Kreuz","und","Becher", "anzunehmen"]).reshape(-1,1)
print(g.test(x))
