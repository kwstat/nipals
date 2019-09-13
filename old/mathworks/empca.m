function [u, s, v, a] = empca(a, ncomps, emtol, maxiters)
%EMPCA	Expectation-Maximization Principal Component Analysis
%   [U, S, V] = EMPCA(A,N) calculates N principal components of matrix A,
%   and returns U, S, V that approximate the N-rank truncation of the
%   singular value decomposition of A. S is a diagonal matrix with singular
%   values corresponding to the contribution of each principal component to
%   matrix A. U and V have orthonormal columns. Matrix A is interpreted as
%   a 2D array. A can be a full or sparse matrix of class 'single',
%   'double', or 'gpuArray'. N must be a positive integer, and is reduced
%   to the minimum dimension of A if higher.
% 
%   [U, S, V, E] = EMPCA(A,N) also returns the residual error matrix E
%   resulting from the PCA decomposition, such that A == U*S*V' + E.
% 
%   [...] = EMPCA(A,N,TOL) keeps principal components when changes in U
%   during the last EM iteration are smaller than TOL, instead of the
%   default value of 1e-6. TOL must be a scalar.
% 
%   [...] = EMPCA(A,N,TOL,MAXITER) keeps principal components after MAXITER
%   EM iterations if convergence has not been reached. If omitted, a
%   maximum of 100 EM iterations are computed. MAXITER must be a positive
%   integer.
% 
%   This function implements the expectation maximization principal
%   component analysis algorithm by Stephen Bailey, available in 
%   http://arxiv.org/pdf/1208.4122v2.pdf
% 
%   Bailey, Stephen. "Principal Component Analysis with Noisy and/or
%   Missing Data." Publications of the Astronomical Society of the Pacific
%   124.919 (2012): 1015-1023.

% keep notation from a = u*s*v'.
% regarding empca paper notation:
%   x   : a
%   phi : u
%   c   : sv
    
a = [50 67 90 98 120;
     55 71 93 102 129;
     65 76 95 105 134;
     50 80 102 130 138;
     60 82 97 135 151;
     65 89 106 137 153;
     75 95 117 133 155]
     
ncomps = 4
%% parameters
if ~exist('emtol','var')
    emtol = 1e-6; % an eigenvector is found when max change in absolute eigenvector difference between EM iterations is below emtol
end
if ~exist('maxiters','var')
    maxiters = 100; % or when maxiters is reached, whatever first
end


emtol = max(emtol,eps(class(a)));  % make sure emtol is not below eps

a = reshape(a,size(a,1),[]); % force it to be 2D
ncomps = min(ncomps,min(size(a))); % reduce number of components if higher than maximum possible rank

u  = @zeros(size(a,1),ncomps,class(a)); % allocate memory for results
sv = @zeros(size(a,2),ncomps,class(a));

normc = @(m)bsxfun(@rdivide,m,sqrt(sum(m.^2))); % returns normalized columns

%% empca
for comp = 1:ncomps
    u(:,comp) = normc(@randn([size(a,1) 1],class(a)));
    for iter = 1:maxiters % repeat until u does not change or too many iterations
        u0 = u(:,comp); % store last iteration's u for comparison
        sv(:,comp) = a'*u(:,comp); % E-step
        u(:,comp) = normc(a*sv(:,comp)); % M-step
        if max(abs(u0-u(:,comp))) <= emtol % check convergence
            break % iter
        end
    end
    disp(['comp ' num2str(comp) ' kept after ' num2str(iter) ' iterations'])
    a = a - u(:,comp)*sv(:,comp)'; % update a removing converged component and leaving residual
end
s = diag(sqrt(sum(sv.^2)));
v = normc(sv);