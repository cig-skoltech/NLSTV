This package includes Matlab scripts that implement the proximal operator of the Non-Local Structure Tensor Total Variation (NL-STV) functional described in the paper:

S. Lefkimmiatis and S. Osher, “Non-local Structure Tensor Functionals for Image Regularization”, IEEE Trans. Computational Imaging, vol. 1, issue 1, pp. 16-29, June 2015

The denoise_demo.m file includes two examples of the proximal map of NL-STV used for denoising a grayscale and a color image, respectively. 

The deconv_demo.m file includes two examples of image deblurring a grayscale and a color image under NL-STV regularization. 


The main routines prox_STVNL_AL and deconv_STVNL_AL depend on the mex script proxSpMat2xNc.c(cpp). In the folder ./mexFunctions/source there are precompiled version for Mac OS and Linux. If you wish to run the script proxSpMat2xNc in a different operating system you will have first to compile the .c/.cpp files. One possible option to do this is by executing the compile.m file located in the ./mexFunctios/source subfolder. 
