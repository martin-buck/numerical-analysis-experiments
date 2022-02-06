% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 2/4/22
%

% construct a matrix that is singular and compute its LU decomposition
A = [.1 .2 .3; .4 .5 .6; .7 .8 .9];
[L, U] = lu(A);

% now solve a system using Gaussian Elimination: Ly=b and Ux=y
b = [.1; .3; .5];
y = L\b;
x = U\y;
