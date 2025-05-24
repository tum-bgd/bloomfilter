# BloomFilter - A Simple and Fast Bloom Filter for Numpy and Python

This repository contains a Bloom Filter implementation for Python. In contrast to existing implementations, we
optimized this one for being fast and flexible for many big data situations.

Therefore, we delibaretly do not work with Python objects but rely on numpy to manage collections of objects. In this way, the number of API calls can be greatly reduced. Furthermore, we do not use the internal hashing of Python but rely on the fast and high-quality Murmur3 hash. 



# Building and Installing
```
$ python3 setup.py bdist_wheel
$ pip3 install dist/*
```

# Quick Docs

- A BF is a class BloomFilter
- core functionality first
  - configure(k, m) sets k and m
  - insert (M, axis=0) inserts each row of M (along axis)
  - test (M, axis=0) tests each row of M
  - clear does what is shall do
  - tobytes: returns a compact bytes repr of the filter contents excluding config
  - frombytes: reads the repr
   
## Data Types for Numpy
The approach we take does work well, however, has a certain risk when not all bytes reserved for objects in
the numpy array are reproducibly initialized. We expect and checked the following:
- for native numeric types, the size reserved is the actual size and, thus, the hashed data is the data
- for strings, when using a char table (e.g., 'U<6'), it works well, because numpy is initializing unused chars to zero

There are still some caveats:
- Bloom filters do hashing on a bit-level, thus, are not compatible when moving BFs between little endian and big endian systems

For whatever you do other than above and to cross-check, you can use the `dbg_getrow` function to get a representation of the actually hashed bytes. Make sure, they are what you expect.
