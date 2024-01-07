N = 48;  % Size of A
% Case 1 
A1 = LSexample(N, 1, 5, 1.0, false);      
b1 = ones(N,1);                           
b1 = b1/sqrt(N);                          
v1 = diag(A1);
sol1 = b1./v1;
% Case 2
A2 = LSexample(N, 1, 100, 1.0, false);
b2 = ones(N,1);
b2 = b2/sqrt(N);
v2 = diag(A2);
sol2 = b2./v2;
% Case 3
A3 = LSexample(N, 1, 5, 0.1, false);
b3 = ones(N,1);
b3 = b3/sqrt(N);
v3 = diag(A3);
sol3 = b3./v3; 
% Case 4
A4 = LSexample(N, 1, 5, 1.0, false);
t = 10^(-13)*ones(N-2,1);
b4 = [1/sqrt(2);t;1/sqrt(2)];
v4 = diag(A4);
sol4 = b4./v4;

tol = 10^(-10);

% relative error and K-bound for Case 1 
[x1, flag1, relres1, iter1, resvec1, ritz1, xvec1] = ConjugateGradient(A1, b1, tol);
[l1, L11] = size(xvec1);
y11 = [];
y21 = [];
t1 = 0;
t2 = 0;
for i = 1:L11
    t1 = Anorm(A1, sol1 - xvec1(:, i)) / Anorm(A1, sol1);
    y11 = [y11, t1];
    t2 = Kbound(A1, i-1);
    y21 = [y21, t2];
end

% relative error and K-bound for Case 2
[x2, flag2, relres2, iter2, resvec2, ritz2, xvec2] = ConjugateGradient(A2, b2, tol);
[l1, L12] = size(xvec2);
y12 = [];
y22 = [];
t1 = 0;
t2 = 0;
for i = 1:L12
    t1 = Anorm(A2, sol2 - xvec2(:, i)) / Anorm(A2, sol2);
    y12 = [y12, t1];
    t2 = Kbound(A2, i-1);
    y22 = [y22, t2];
end

% relative error and K-bound for Case 3 
[x3, flag3, relres3, iter3, resvec3, ritz3, xvec3] = ConjugateGradient(A3, b3, tol);
[l1, L13] = size(xvec3);
y13 = [];
y23 = [];
t1 = 0;
t2 = 0;
for i = 1:L13
    t1 = Anorm(A3, sol3 - xvec3(:, i)) / Anorm(A3, sol3);
    y13 = [y13, t1];
    t2 = Kbound(A3, i-1);
    y23 = [y23, t2];
end

% relative error and K-bound for Case 4
[x4, flag4, relres4, iter4, resvec4, ritz4, xvec4] = ConjugateGradient(A4, b4, tol);
[l1, L14] = size(xvec4);
y14 = [];
y24 = [];
t1 = 0;
t2 = 0;
for i = 1:L14
    t1 = Anorm(A4, sol4 - xvec4(:, i)) / Anorm(A4, sol4);
    y14 = [y14, t1];
end
for i = 1:20
    t2 = Kbound(A4, i-1);
    y24 = [y24, t2];
end

figure

subplot(2,2,1)
semilogy(y11, 'b--')
hold on
semilogy(y21, 'k-.')
hold on
vlastc1 = eig(A1);
for  l = 2 : 2 : L11 + 1           % computing of bound based on removing of outlying eigenvalues, Case 1
   LiStb1 = [];
   kapal = vlastc1(N - l + 1)/vlastc1(1); 
   for m = l : L11 
        temp1 = 2 * ((sqrt(kapal) - 1)/(sqrt(kapal) + 1))^(m - l); 
        LiStb1 = [LiStb1,temp1];
   end
 semilogy(l:L11,LiStb1);
 hold on
end
title('Case 1');
hold off

subplot(2,2,2)
semilogy(y12, 'b--')
hold on
semilogy(y22, 'k-.')
hold on
vlastc = eig(A2);
for l = 2 : 2 : L12 + 1          % computing of bound based on removing of outlying eigenvalues, Case 2
   LiStb = [];
   kapal = vlastc(N - l + 1)/vlastc(1); 
   for m = l : L12  
        temp1 = 2 * ((sqrt(kapal) - 1)/(sqrt(kapal) + 1))^(m - l); 
        LiStb = [LiStb,temp1];
   end
   semilogy(l:L12,LiStb);
   hold on
end
title('Case 2'); 
hold off

subplot(2,2,3)
semilogy(y13,'b--');
hold on
semilogy(y23,'k-.');
vlastc = eig(A3);
for l = 2 : L13 + 1       % computing of bound based on removing of outlying eigenvalues, Case 3
   LiStb = [];
   kapal = vlastc(N - l + 1)/vlastc(1);  
   for m = l : L13        
        temp1 = 2 * ((sqrt(kapal) - 1) / (sqrt(kapal) + 1))^(m - l); 
        LiStb = [LiStb, temp1];
   end
   semilogy(l:L13,LiStb);
   hold on
end
title('Case 3');
hold off

subplot(2,2,4)
semilogy(y14,'b--')
hold on
semilogy(y24, 'k-.')
hold on
vlastc = eig(A4);
for l = 2 : L14 + 1          % computing of bound based on removing of outlying eigenvalues, Case 4
   LiStb = [];
   kapal = vlastc(N - l + 1)/vlastc(1); 
   for m = l : 20
        temp1 = 2 * ((sqrt(kapal) - 1)/(sqrt(kapal) + 1))^(m - l);      
        LiStb = [LiStb,temp1];
   end
   semilogy(l:20,LiStb);
   hold on
end
title('Case 4');
hold off