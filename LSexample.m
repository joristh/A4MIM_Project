function [A, b] = LSexample(N, l1, lN, gam, reverse)
% modified Liesen Strakos (2013) example
% Input
% N : size of A. Default N = 34
% l1,lN : smallest and largest eigenvalue. Default l1=0.1, lN=100
% gam :   gam < 1 most of the eigenvalues on the left-hand side of the spectrum
%         gam = 1 equispaced eigenvalues
%         Default gam=0.65
% reverse : true/false to accumulate eigenvalues towards the largest
%
% Output
% A : square diagonal matrix
% b : random right-hand side

if nargin  == 0
    N = 34;
    l1 = 0.1;
    lN = 100;
    gam = 0.65;
    reverse = false;
end

if nargin  == 1
    l1 = 0.1;
    lN = 100;
    gam = 0.65;
    reverse = false;
end

if nargin == 3
    gam = 0.65;
    reverse = false;
end

if nargin == 4
    reverse = false;
end

if reverse == false
    l(1) = l1;
    for j=2:N-1
        l(j) = l1 + (j-1)/(N-1)*(lN-l1)*gam^(N-j);
    end
    l(N) = lN;
else
    l(1) = lN;
    for j=2:N-1
        l(j) = lN - (j-1)/(N-1)*(lN-l1)*gam^(N-j);
    end
    l(N) = l1;
end

A = diag(l);
b = rand(N,1);
