# Matrix Multiplication

> performance optimization via register use and cache reuse in dynamically allocated matrix multiplication

## Instruction
- `make` compile all source files
- `make clean` remove outputted files
- `make test` execute the cli
- `make valgrind` run valgrind

## Algorithms

### Single Register Use
```C
for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
        register double r = c[i*n+j];
        for (k=0; k<n; k++)
            r += a[i*n+k] * b[k*n+j];
        c[i*n+j];
    }
```

### Blocked Version Algorithm
```C
for (i = 0; i < n; i+=B)
    for (j = 0; j < n; j+=B)
        for (k = 0; k < n; k+=B)
            for (i1 = i; i1 < i+B; i1++)
                for (j1 = j; j1 < j+B; j1++) {
                    register double r=c[i1*n+j1];
                        for (k1 = k; k1 < k+B; k1++)
                            r += a[i1*n + k1]*b[k1*n + j1]; 
                        c[i1*n+j1]=r; 
                    }
```

## Course Information
<pre align="center">
CS 481 - High Performance Computing
Instructor: Dingwen Tao
</pre>