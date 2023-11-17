clear all
N = 30;
A = LSexample(N, 0.1, 10e3, 0.6);
b = ones(N, 1)/sqrt(N);

figure(1)
clf
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

semilogy(1:iter, err_norm/xnorm)
hold on
semilogy(1:iter, err_norm_ro/xnorm)
ylim([1e-8, 10])
legend("CG", "CGRO")