function compile_mex 

if ispc
  mex -v proxSpMat2xNc.cpp libmwblas.lib CFLAGS="\$CFLAGS /openmp" LDFLAGS="\$LDFLAGS /openmp" -largeArrayDims  -I../headers/
else
  mex -v proxSpMat2xNc.c -lmwblas CFLAGS="\$CFLAGS -fopenmp -std=c99“ LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims  -I../headers/
end
