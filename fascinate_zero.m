function [ f ] = fascinate_zero( A_row,F,alpha,beta)
%COLDSTARTINC Summary of this function goes here
%   coldStartInc function incrementally update teh low-rank matrix of layer
%   i in the mulan. The update rule is as F_new(n+1,t) =
%   {alphai*A'(1:n,n+1)*F(:,t)}/{beta+alphai*(deg(n+1)-A(n+1,n+1))}
n = size(F,1);
numerators = alpha*A_row(1:n)*F;
deg_new = sum(A_row);
denominators = beta + alpha*(deg_new-A_row(n+1));
f = numerators/denominators;


end

