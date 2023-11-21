function [x, flag, relres, iter, resvec, ritz] = ConjugateGradient(A, b, tol, maxit, x0)
%
%   ConjugateGradient tries to solve the linear system A*x = b for x iteratively with the conjugate gradients method where A is a symmetric positive definite matrix.
%   When the attempt is successful, Conjugate Gradient displays a message to confirm convergence.
%   If ConjugateGradient fails to converge after the maximum number of iterations or halts for any reason,
%   it displays a diagnostic message that includes the relative residual norm(b-A*x)/norm(b) and the iteration number at which the method stopped.
%   Displayed messages are omitted if the flag is returned
%   
%   POSSIBLE FUNCTION SIGNATURES
%
%   x = ConjugateGradient(A,b)
%   x = ConjugateGradient(A,b,tol)
%   x = ConjugateGradient(A,b,tol,maxit)
%   x = ConjugateGradient(A,b,tol,maxit,x0)
%   [x,flag] = ConjugateGradient(___)
%   [x,flag,relres] = ConjugateGradient(___)
%   [x,flag,relres,iter] = ConjugateGradient(___)
%   [x,flag,relres,iter,resvec] = ConjugateGradient(___)
%
%   INPUT ARGUMENTS
%   
%   A      - symmetric positive definite n x n matrix ()
%   b      - right hand side 1 x n vector
%   tol    - tolerance on relative residual norm(r)/norm(b) (default 1e-6)
%   maxit  - maximum number of iterations (default min(n, 20))
%   x0     - initial guess (default zero vector)
%   
%   OUTPUT ARGUMENTS
%   
%   x      - approximate solution
%   flag   - specifies whether the algorithm successfully converged. When flag = 0, convergence was successful
%   relres - relative residual norm(b-A*x)/norm(b). If flag is 0, then relres <= tol
%   iter   - iteration number iter at which x was computed
%   resvec - vector of the residual norm at each iteration, including the first residual norm(b-A*x0)
%   ritz   - vector of Ritz values
%
%   FLAG VALUES
%   
%   0      - Success: ConjugateGradient converged to the desired tolerance tol within maxit iterations.
%   1      - Failure: ConjugateGradient iterated maxit iterations but did not converge.
%   4      - Failure: One of the scalar quantities calculated by the ConjugateGradient algorithm became too small or too large to continue computing.
%   

% check input arguments
try
    % default parameters
    if nargin < 2
        errID = 'ConjugateGradient:inputError';
        msgtext = 'Not enough inputs. Needs at least matrix and rhs.';
        throw(MException(errID, msgtext));
    end
    if nargin < 3
        tol = 1e-6;
    end
    if nargin < 4
        maxit = min(size(A,1), 20);
    end
    if nargin < 5
        x0 = zeros(length(b), 1);
    end
    
    % check if matrix is square
    [n, m] = size(A);
    if n ~= m
        errID = 'ConjugateGradient:inputError';
        msgtext = 'Input matrix is not square. Size is %d x %d.';
        throw(MException(errID, msgtext, n, m));
    end
    
    % check if matrix is symmetric
    if ~all(all(A == A'))
        errID = 'ConjugateGradient:inputError';
        msgtext = 'Input matrix is not symmetric';
        throw(MException(errID, msgtext));
    end
    
    % check rhs
    if length(b) ~= n
        errID = 'ConjugateGradient:inputError';
        msgtext = 'Length of rhs %d does not match matrix size %d.';
        throw(MException(errID, msgtext, length(b), n));
    end
    
    % check initial guess
    if length(x0) ~= n
        errID = 'ConjugateGradient:inputError';
        msgtext = 'Length of initial guess %d does not match matrix size %d.';
        throw(MException(errID, msgtext, length(x0), n));
    end
catch exception
    disp('ConjugateGradient:inputError')
    throw(exception)
end

alpha = 0;
beta = 0;
r = b - A*x0;
r2 = r' * r;
p = r;
Ap = zeros(1, n);
x = x0;
resvec_ = [r2];

flag_ = -1;
abstol2 = (tol*norm(b))^2;

if nargout > 5
    residuals = r/sqrt(r2);
    store_residuals = true;
else
    store_residuals = false;
end

try
    for i=1:maxit
        iter_ = i;
        Ap = A*p;
        alpha = r2 / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        r2 = r'*r;
        resvec_(i+1) = r2;
        if store_residuals
            % storing orthonormal basis for Ritz values
            residuals = [residuals, r/sqrt(r2)];
        end
        if r2 <= abstol2
            % success: method converged to specified tolerance
            flag_ = 0;
            relres_ = norm(r)/norm(b);
            if nargout <= 1
                disp(['ConjugateGradient converged at iteration ', num2str(iter_), ' to a solution with relative residual ', num2str(relres_), '.'])
            end
            break
        end
        beta = r2/resvec_(i);
        p = r + beta * p;
    end
catch
    % failure: numerical exception during iteration
    flag_ = 4;
    iter_ = 0;
    relres_ = 0;
    if nargout <= 1
        disp('One of the scalar quantities calculated by the ConjugateGradient algorithm became too small or too large to continue computing.')
    end
end

if flag_ < 0
    % failure: maximum number of iterations reached
    flag_ = 1;
    iter_ = maxit;
    relres_ = norm(r)/norm(b);
    if nargout <= 1
        disp(['ConjugateGradient stopped at iteration ', num2str(maxit), ' without converging to the desired tolerance ', num2str(tol)])
        disp('because the maximum number of iterations was reached.')
        disp(['The iterate returned (number ', num2str(maxit), ') has relative residual ', num2str(relres_),'.'])
    end
end

% variable number of output arguments
if nargout > 1
    flag = flag_;
end
if nargout > 2
    relres = relres_;
end
if nargout > 3
    iter = iter_;
end
if nargout > 4
    resvec = sqrt(resvec_(:));
end
if nargout > 5
    ritz = eig(residuals'*A*residuals);
end