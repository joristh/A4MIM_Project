clear all

% set up problem
N = 30;
A = LSexample(N, 0.1, 1e3, 0.6);
b = ones(N, 1)/sqrt(N);

% compare convergence in exact and finite arithmetic 
iter = 100;
err_norm = zeros(iter, 1);
err_norm_ro = zeros(iter, 1);
x = A\b;
xnorm = x'*A*x;

for i = 1:iter
    [xcg, flag] = ConjugateGradient(A, b, 1e-6, i);
    [xcgro, flag] = ConjugateGradientRO(A, b, 1e-6, i);
    err1 = x - xcg;
    err2 = x - xcgro;
    err_norm(i) = err1'*A*err1;
    err_norm_ro(i) = err2'*A*err2;
end

figure(1)
clf
semilogy(1:iter, err_norm/xnorm)
hold on
semilogy(1:iter, err_norm_ro/xnorm)
ylim([1e-8, 10])
legend("CG (finite arithmetic)", "CGRO (exact arithmetic)")
title("Convergence comparison in A-norm")

% cumulative spectral density plots of eigenvalues and Ritz values
[xcgro, flag, relres, iter, resvec, ritz10] = ConjugateGradientRO(A, b, 1e-6, 10);
[xcgro, flag, relres, iter, resvec, ritz20] = ConjugateGradientRO(A, b, 1e-6, 20);

figure(2)
clf
spectraldensity(diag(A))
hold on
spectraldensity(ritz10)
spectraldensity(ritz20)
xlim([0, 1e3])
legend("eigenvalues", "Ritz values k=10", "Ritz values k=20")
legend('Location','southeast')
title("Eigenvalues and Ritz values (exact)")
