function [ D_infer ] = restoreD( G_new, F )
%EVALAUC Summary of this function goes here
%   Detailed explanation goes here
uniqueMap = triu(G_new);
[I,J] = find(uniqueMap);

for i = 1:length(I)
    index = uniqueMap(I(i),J(i)); 
    U = F{I(i)}.F;
    V = F{J(i)}.F;
    D_infer{index}.D = U*V';
end

end

