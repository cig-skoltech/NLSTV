%Function proxSpMat2xNc solves the following matrix regularization problem 
%
% min 1/2 ||X-Y||^2 + rho ||X||_p
%  X
%
% where Y is a 2xNc real matrix and p is the order of the corresponding
% Schatten norm.
%
%Matlab Usage: X=proxSpMat2xNc(Y,p,rho,c0);
% c0:initial guess for the root.
%
% Function proxSpMat2xNc can handle multiple matrix inputs.
% Y: Y is a multidimensional matrix where its two last dimensions should be 
% of size 2xNc.
%
% Note that in the multiple matrix input case we can also use a
% different rho for each matrix. Then numel(rho) should be either equal to
% the number of the input matrices or should be a scalar, which means
% that all the matrices are using the same rho.
% =========================================================================
%
%  Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================