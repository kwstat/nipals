% use this file

function [u, s, v, a] = empca_w(a, w, ncomp, emtol, maxiters)
%EMPCA	Expectation-Maximization Principal Component Analysis
%   [U, S, V] = EMPCA(A,W,N) calculates N principal components of matrix A,
%   using weight matrix W.
%   Returns U, S, V that approximate the N-rank truncation of the singular
%   value decomposition of A. S is a diagonal matrix with singular values
%   corresponding to the contribution of each principal component to matrix
%   A. U and V have orthonormal columns. Matrix A is interpreted as a 2D
%   array. A can be a full or sparse matrix of class 'single', 'double', or
%   'gpuArray'. N must be a positive integer, and is reduced to the minimum
%   dimension of A if higher.
% 
%   [U, S, V, E] = EMPCA(A,W,N) also returns the residual error matrix E
%   resulting from the PCA decomposition, such that A == U*S*V' + E.
% 
%   [...] = EMPCA(A,W) calculates all the components.
% 
%   [...] = EMPCA(A,W,N,TOL) keeps principal components when changes in U
%   during the last EM iteration are smaller than TOL, instead of the
%   default value of 1e-6. TOL must be a scalar.
% 
%   [...] = EMPCA(A,W,N,TOL,MAXITER) keeps principal components after MAXITER
%   EM iterations if convergence has not been reached. If omitted, a
%   maximum of 100 EM iterations are computed. MAXITER must be a positive
%   integer.
% 
%   This function implements the expectation maximization principal
%   component analysis algorithm by Stephen Bailey, available in 
%   http://arxiv.org/pdf/1208.4122v2.pdf
%   Bailey, Stephen. "Principal Component Analysis with Noisy and/or
%   Missing Data." Publications of the Astronomical Society of the Pacific
%   124.919 (2012): 1015-1023.
% regarding empca paper notation:
%   phi : u
%   c   : C

% ------------------

B1 = [50 67 90 98 120;
      55 71 93 102 129;
      65 76 95 105 134;
      50 80 102 130 138;
      60 82 97 135 151;
      65 89 106 137 153;
      75 95 117 133 155];
B1wt = [1 1 1 1 1;
        1 1 1 1 1;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ];
% In B2, the first two elements of the first column are really NA,
% but the code will set the values to 0 and then use 0 weight
B2 = [0 67 90 98 120;
      0 71 93 102 129;
      65 76 95 105 134;
      50 80 102 130 138;
      60 82 97 135 151;
      65 89 106 137 153;
      75 95 117 133 155];
B2wt = [0 1 1 1 1;
        0 1 1 1 1;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ;
        1 1 1 1 1 ];

% -----------------------------
     
X = B2
W = B2wt
ncomp = 4;


% column-centered & scaled     
X = bsxfun(@minus, X, mean(X));
X = bsxfun(@rdivide, X, std(X))
 
%% parameters to determine when an eigenvector is found: 
emtol = 1e-12; % max change in absolute eigenvector difference between EM iterations
maxiters = 100; % or when maxiters is reached, whatever first

emtol = max(emtol,eps(class(X)));  % make sure emtol is not below eps

%X = reshape(X,size(X,1),[]); % force it to be 2D

if ~exist('ncomp','var')
    ncomp = min(size(X)); % set to max rank if not specified
else
    warning 'ncomp reduced to max rank'
    ncomp = min(ncomp,min(size(X))); % reduce if higher than maximum possible rank
end

% allocate space for results
P  = @zeros(size(X,1),ncomp,class(X));
C = @zeros(size(X,2),ncomp,class(X));

% macrp tp return normalized columns
normc = @(m)bsxfun(@rdivide,m,sqrt(sum(m.^2))); 

%W = ~isnan(X);
%X(~W) = 0;

%% empca
for comp = 1:ncomp
    % random direction vector
    P(:,comp) = normc(@randn([size(X,1) 1],class(X)));
    for iter = 1:maxiters % repeat until u does not change
        P0 = P(:,comp); % store previous u for comparison
        
        % E step
        C(:,comp) = X' * P(:,comp);

        % M-step with weights. Bailey eqn 21. sum(,2)= rowsums
        % .* is element-wise
        % column C[,h] times each column of t(W)
        CW = bsxfun(@times, C(:,comp) , W');
        P(:,comp) = sum(X .* CW',2) ./ (CW' * C(:,comp));
        P(:,comp) = normc(P(:,comp));
        
        if max(abs(P0-P(:,comp))) <= emtol
            break
        end
    end
    disp(['comp ' num2str(comp) ' kept after ' num2str(iter) ' iterations'])

    X = X - P(:,comp) * C(:,comp)'; % deflate X = X - P*C'
    X(~W) = 0;
    disp("test")
end

% restore missing values
% X(~W) = NaN;
s2 = diag(sum(C.^2));
v = normc(C);

disp("-----------------------------");
output_precision(2)
disp ("Value of P (scores)"), disp (eval(mat2str(P,3)));
disp("s2"), disp(s2);
disp("v (loadings)"), disp(v);
disp("Value of C"), disp(C);
disp("P'P"), disp(round(P'*P)) % check U is orthonormal
