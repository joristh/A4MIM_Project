clear all

% set up
n = 10;
N = 100;
cluster_spacing = 1e-12;
cl_mod = (1-n/2 : n/2) * cluster_spacing;
b = ones(N, 1)/sqrt(N);

A1 = LSexampleMod(n, .1, 1e3, .6, true);
A2 = LSexampleMod(n, .1, 1e3, .6, false);
A3 = LSexampleMod(n, .1, 1e3, 1, false);
a1 = diag(A1);
a2 = diag(A2);
a3 = diag(A3);
lambda1 = [];
lambda2 = [];
lambda3 = [];

for i = 1:n
    lambda1 = [lambda1, a1(i) + cl_mod];
    lambda2 = [lambda2, a2(i) + cl_mod];
    lambda3 = [lambda3, a3(i) + cl_mod];
end

A1 = diag(lambda1);
A2 = diag(lambda2);
A3 = diag(lambda3);

%%
% solution
x1 = A1\b;
x2 = A2\b;
x3 = A3\b;

% errors for exact CG
x1norm = sqrt(x1'*A1*x1);
x2norm = sqrt(x2'*A2*x2);
x3norm = sqrt(x3'*A3*x3);

iter = 20;
err1 = zeros(iter, 1);
err2 = zeros(iter, 1);
err3 = zeros(iter, 1);
err1(1) = x1norm;
err2(1) = x2norm;
err3(1) = x3norm;

for i = 1:iter
    [X1, ~, ~, ~, ~, ritz1] = ConjugateGradientRO(A1, b, 1e-10, i);
    [X2, ~, ~, ~, ~, ritz2] = ConjugateGradientRO(A2, b, 1e-10, i);
    [X3, ~, ~, ~, ~, ritz3] = ConjugateGradientRO(A3, b, 1e-10, i);
    e1 = x1 - X1;
    e2 = x2 - X2;
    e3 = x3 - X3;
    err1(i+1) = sqrt(e1'*A1*e1);
    err2(i+1) = sqrt(e2'*A2*e2);
    err3(i+1) = sqrt(e3'*A3*e3);

    % ritz values at certain iterations
    if i==4
        ritz5_1 = ritz1;
        ritz5_2 = ritz2;
        ritz5_3 = ritz3;
    end
    if i==9
        ritz10_1 = ritz1;
        ritz10_2 = ritz2;
        ritz10_3 = ritz3;
    end
    if i==14
        ritz15_2 = ritz2;
    end
end

% plot of relative errors
figure(1)
semilogy(0:20, err1/x1norm, "b");
hold on
semilogy(0:20, err2/x2norm, "r");
hold on
semilogy(0:20, err3/x3norm, "g");
ylim([1e-11,10])
hold off

%%
% cumulative spectral density
[csd_l1, l1] = ecdf(lambda1);
[csd_l2, l2] = ecdf(lambda2);
[csd_l3, l3] = ecdf(lambda3);

[csd_ritz5_1, ritz5_1] = ecdf(ritz5_1);
[csd_ritz10_1, ritz10_1] = ecdf(ritz10_1);

[csd_ritz5_2, ritz5_2] = ecdf(ritz5_2);
[csd_ritz10_2, ritz10_2] = ecdf(ritz10_2);
[csd_ritz15_2, ritz15_2] = ecdf(ritz15_2);

[csd_ritz5_3, ritz5_3] = ecdf(ritz5_3);
[csd_ritz10_3, ritz10_3] = ecdf(ritz10_3);

% CSD plots
figure(2)
subplot(3,1,1)
stairs(l1,csd_l1,"Color",[0, .7, 1],"LineWidth",2,"DisplayName","eigevalues")
hold on
stairs(ritz5_1,csd_ritz5_1,":","Color",[0,0,0],"DisplayName","Ritz values for k=5")
hold on
stairs(ritz10_1,csd_ritz10_1,"--","Color",[0,0,0],"DisplayName","Ritz values for k=10")
hold off
legend("Location","northwest")
title('accumulated to the right')
xlim([0,1001])

subplot(3,1,2)
stairs(l2,csd_l2,"Color",[1, .35, .35],"LineWidth",2,"DisplayName","eigenvalues")
hold on
stairs(ritz5_2,csd_ritz5_2,":","Color",[0,0,0],"DisplayName","Ritz values for k=5")
hold on
stairs(ritz10_2,csd_ritz10_2,"-.","Color",[0,0,0],"DisplayName","Ritz values for k=10")
hold on
stairs(ritz15_2,csd_ritz15_2,"--","Color",[0,0,0],"DisplayName","Ritz values for k=15")
hold off
legend("Location","southeast")
title('accumulated to the left')
xlim([0,1001])

subplot(3,1,3)
stairs(l3,csd_l3,"g","LineWidth",2,"DisplayName","eigenvalues")
hold on
stairs(ritz5_3,csd_ritz5_3,":","Color",[0,0,0],"DisplayName","Ritz values for k=5")
hold on
stairs(ritz10_3,csd_ritz10_3,"--","Color",[0,0,0],"DisplayName","Ritz values for k=10")
hold off
legend("Location","southeast")
title('equally spaced')
xlim([0,1001])
