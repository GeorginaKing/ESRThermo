function corrR = rankCorr(Cstar,R)
%CORRR = RANKCORR(CSTAR,R)
% Induces rank correlation onto a sample.
% Method of Iman & Conover.
% Iman, R. L., and W. J. Conover. 1982.
% ?A Distribution-free Approach to Inducing Rank
% Correlation Among Input Variables.?
% Communications in Statistics: Simulation and
% Computations 11:311-334.
% Input:
% Cstar     : wanted correlation matrix (k,k)
% R         : sample matrix (N,k)
% Output:
% corrR     : correlated sample matrix (N,k)
%
C = Cstar;
[N k] = size(R);
% Calculate the sample correlation matrix T
T = corrcoef(R);
% Calculate lower triangular cholesky
% decomposition of Cstar
% i.e. P*P? = C
P = chol(C)';
% Calculate lower triangular cholesky decomposition
% of T, i.e. Q*Q? = T
Q = chol(T)';
% S*T*S? = C
S = P*inv(Q);
% Replace values in samples with corresponding
% rank-indices and convert to van der Waerden scores
RvdW = -sqrt(2).*erfcinv(2*(getRanks(R)/(N+1)));
% Matrix RBstar has a correlation matrix exactly
% equal to C
RBstar = RvdW*S';
% Match up the rank pairing in R according to RBstar
ranks = getRanks(RBstar);
sortedR = sort(R);
for i=1:k
    corrR(:,i) = sortedR(ranks(:,i),i);
end

function r = getRanks(u)
%R = GETRANKS(U)
% Ranking of a vector (matrix)
% Input:
%   u : vector (matrix) (nrow,ncol)
% Output:
%   r : rank of the vector (nrow, ncol)
% Returns a matrix the size of u where each
% element in u has been replaced by it?s
% corresponding rank value.
%
% If a matrix is used as argument, each column
% is treated separately.
%
if size(u,1)==1
    u=u';
end
[nr nc] = size(u);
[s,ind] = sort(u);
for i=1:nc
    r(ind(1:nr,i),i) = 1:nr;
end