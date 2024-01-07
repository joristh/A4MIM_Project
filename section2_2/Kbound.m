function [bo] = Kbound(A, k)
K = cond(A);
t1 = sqrt(K);
bo = 2 * ((t1 - 1) / (t1 + 1)) ^ k;