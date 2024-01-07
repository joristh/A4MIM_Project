function [An] = Anorm(A, x)
An = sqrt(x' * A * x);