Past exams code solutions

**Requirements**
```
module load gcc-glibc
module load dealii
```

**Compile and Run**
```
cd build
cmake ..
make
./executable
```

If it gives any Cmake errors:

Inside the build dir
```
rm CMakeCache.txt
```
