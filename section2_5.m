clear all

% finite arithmetic
N = 30;

A1 = LSexample(N, 0.1, 1e3, 0.6, true);
A2 = LSexample(N, 0.1, 1e3, 0.6);
A3 = LSexample(N, 0.1, 1e3, 1.0);

b = ones(N, 1)/sqrt(N);

x1 = A1\b;
x2 = A2\b;
x3 = A3\b;

x1norm = x1'*A1*x1;
x2norm = x2'*A2*x2;
x3norm = x3'*A3*x3;

iter = 50;
err1 = zeros(iter, 1);
err2 = zeros(iter, 1);
err3 = zeros(iter, 1);
err1(1) = x1norm;
err2(1) = x2norm;
err3(1) = x3norm;

for i = 1:iter
    [X1, flag] = ConjugateGradient(A1, b, 1e-10, i);
    [X2, flag] = ConjugateGradient(A2, b, 1e-10, i);
    [X3, flag] = ConjugateGradient(A3, b, 1e-10, i);
    e1 = x1 - X1;
    e2 = x2 - X2;
    e3 = x3 - X3;
    err1(i+1) = e1'*A1*e1;
    err2(i+1) = e2'*A2*e2;
    err3(i+1) = e3'*A3*e3;
end

% error convergence comparison finite precision
figure(1)
clf
semilogy(0:iter, err1/x1norm, "b", linewidth=2)
hold on
semilogy(0:iter, err2/x2norm, "r", linewidth=2)
semilogy(0:iter, err3/x3norm, "g", linewidth=2)
ylim([1e-8, 10])
xlim([0, 50])
exportgraphics(gcf,'plots/sec2_5_finite.pdf','ContentType','vector')

% exact arithmetic modified problem
tmp = diag(A1);
A1_mod = diag([tmp; tmp+1e-14; tmp+2e-14; tmp-1e-14]);
tmp = diag(A2);
A2_mod = diag([tmp; tmp+1e-14; tmp+2e-14; tmp-1e-14]);
tmp = diag(A3);
A3_mod = diag([tmp; tmp+1e-14; tmp+2e-14; tmp-1e-14]);

b_mod = ones(4*N, 1)/sqrt(4*N);

x1 = A1_mod\b_mod;
x2 = A2_mod\b_mod;
x3 = A3_mod\b_mod;

x1norm = x1'*A1_mod*x1;
x2norm = x2'*A2_mod*x2;
x3norm = x3'*A3_mod*x3;

iter = 50;
err1 = zeros(iter, 1);
err2 = zeros(iter, 1);
err3 = zeros(iter, 1);
err1(1) = x1norm;
err2(1) = x2norm;
err3(1) = x3norm;

for i = 1:iter
    [X1, flag] = ConjugateGradient(A1_mod, b_mod, 1e-10, i);
    [X2, flag] = ConjugateGradient(A2_mod, b_mod, 1e-10, i);
    [X3, flag] = ConjugateGradient(A3_mod, b_mod, 1e-10, i);
    e1 = x1 - X1;
    e2 = x2 - X2;
    e3 = x3 - X3;
    err1(i+1) = e1'*A1_mod*e1;
    err2(i+1) = e2'*A2_mod*e2;
    err3(i+1) = e3'*A3_mod*e3;
end

% error convergence comparison finite precision
figure(2)
clf
semilogy(0:iter, err1/x1norm, "b", linewidth=2)
hold on
semilogy(0:iter, err2/x2norm, "r", linewidth=2)
semilogy(0:iter, err3/x3norm, "g", linewidth=2)
ylim([1e-8, 10])
xlim([0, 50])
exportgraphics(gcf,'plots/sec2_5_exact.pdf','ContentType','vector')

%% 
% cumulative spectral density plots of eigenvalues and Ritz values

fig = figure(3)
clf
set(gcf,'units','points','position',[10,10,315,500])

subplot(3,1,1);
stairs(sort(diag(A1)), (1:N)/N, "b", linewidth=1)
hold on
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A1, b, 1e-10, 10);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k:", linewidth=1)
title("acc. to the right")
legend("eigenvalues", "Ritz values k = 10")
legend('Location','northwest')
xlim([0, 1000])
ylim([0, 1])
xticks([0, 200, 400, 600, 800, 1000])
yticks([0, 0.5, 1])

subplot(3,1,2); 
stairs(sort(diag(A2)), (1:N)/N, "r", linewidth=1)
hold on
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A2, b, 1e-10, 10);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k:", linewidth=1)
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A2, b, 1e-10, 20);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k-.", linewidth=1)
title("acc. to the left")
legend("eigenvalues", "Ritz values k = 10", "Ritz values k = 20")
legend('Location','southeast')
xlim([0, 1000])
ylim([0, 1])
xticks([0, 200, 400, 600, 800, 1000])
yticks([0, 0.5, 1])

subplot(3,1,3); 
stairs(sort(diag(A3)), (1:N)/N, "g", linewidth=1)
hold on
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A3, b, 1e-10, 10);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k:", linewidth=1)
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A3, b, 1e-10, 20);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k-.", linewidth=1)
[x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A3, b, 1e-10, 30);
stairs(sort(ritz), (1:length(ritz))/length(ritz), "k--", linewidth=1)
title("equally spaced")
legend("eigenvalues", "Ritz values k = 10", "Ritz values k = 20", "Ritz values k = 30")
legend('Location','northwest')
xlim([0, 1000])
ylim([0, 1])
xticks([0, 200, 400, 600, 800, 1000])
yticks([0, 0.5, 1])

set(findall(fig, 'Type', 'Text'),'FontWeight', 'Normal')

exportgraphics(gcf,'plots/sec2_5_ritzvalues.pdf','ContentType','vector')
