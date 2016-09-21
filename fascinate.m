function [ F ] = fascinate( G,D_map,DM,alpha,beta,weight, rank, paras )
%MULANIMPUTE Summary of this function goes here
%   the function takes a multi-layered network as input and output the low
%   rank approximation for each layer
% \sum||W.(D(i,j)-Fi'Fj))||^2 + \alpha\sum tr(Fi'(DAi-Ai)Fi)+ \beta\sum||Fi||^2
% INPUT: G: adjacency matrix for each layers; D_map: dependency map matrix;
% DM: observed dependency matrice; alpha:coefficient of laplacian regularization; beta:coefficient
% of ||Fi||;
% OUTPUT: set of low-rank approximations

if nargin<8
    % maxIter, threshold
    paras = [100,1e-8];
end

[F] = updateF(G,D_map,DM,alpha,beta,weight, rank,paras);

end

function [F] = updateF(G,D_map,DM,alpha,beta,weight, rank, paras)
%[I,J,V] = find(D_map);
maxIter = paras(1);
threshold = paras(2);

%initialized low rank approximations
for i = 1:length(G)
    m = size(G{i}.A,1);
    F{i}.F =  (rand(m, rank))/sqrt(rank);
end


maxChange = 1;
iter = 0;

%iteratively update each low-rank matrix
while maxChange > threshold && iter < maxIter
    %update each low-rank approximation matrix
    for i = 1:length(G)     
        F0 = F{i}.F;
        F{i}.F = updateFi(i,D_map,weight,G,DM,F,alpha,beta);
        changes(i) = max(max(abs(F{i}.F-F0)));
    end
    maxChange  = max(changes);
    iter = iter + 1;
    fprintf('max element changed: %f, iteration: %d\n',maxChange,iter);
end

end

function [Fi_new] = updateFi(curLayer,Dmap,weight,G,DM,F,alpha,beta)

Fi = F{curLayer}.F;
[m,n] = size(Fi);
%initialize A2 and B2
A2 = zeros(m,n);
B2 = zeros(m,n);
%find related layers
relLayers = find(Dmap(curLayer,:));
%get summation part of A2 and B2
for i = 1:length(relLayers)
    j = relLayers(i);
    index = Dmap(curLayer,j);
    Dij = DM{index}.D;
    if curLayer > relLayers(i)
        Dij = Dij';
    
    end
    A2 = A2 + Dij*F{j}.F;
    
    Rhat = get_UVT(Dij,Fi,F{j}.F);
    B2 = B2 + (1-weight)*Rhat*F{j}.F+weight*Fi*(F{j}.F'*F{j}.F);
end
%remaining part of A2
A2 = A2 + beta*G{curLayer}.A*Fi;
%remaining parts of B2
sumA = sum(G{curLayer}.A,2);
DA = spdiags(sumA,0,m,m);
B2 = B2 + alpha*Fi + beta*DA*Fi + 1e-8;

Fi_new = Fi.*sqrt(A2./B2);

end


function [UVT] = get_UVT(R, U, V)

[m, n] = size(R);

[I, J] = find(R);
iSize = size(I, 1);
K = ones(iSize,1);
count = 0;
for j = 1:n
    id = find(R(:,j));
    len = length(id);
    K(count+1:count+len) = U(id, :) * V(j,:)';
    count = count + len;
end

UVT = sparse(I, J, K, m, n);

end