# Exercise 00

## 03
### Build and Run
1. Create an out of source directory  
    ``mkdir build``
2. Switch into the directory  
    ``cd build``
3. Invoke cmake  
    ``cmake ..``  
    (optionally with flags)  
    ``cmake -DCMAKE_CXX_FLAGS="-O3 -march=native" ..``
4. Build the project  
    ``make``
5. Finally, run it  
    ``/main``


## 04
### Build and Run
1. Build:   
   ```g++ -o daxpy.out daxpy.cc -pthread -std=c++17```
2. Run:   
   ```./daxpy.out```