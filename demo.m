clear
load AminerData
%dataset description
%network is constructed from Aminer dataset 
%G contains all adjacency matrices for each layer
%G{1}.A: author-author collaboration network; 
%G{2}.A paper-paper citation network
%G{3}.A: venue-venue citation network
%G_new: the indexed map of cross-layer dependencies
%DO: observed cross-layer dependency matrices, 50% of the entire data DU
%DU: complete cross-layer dependency matrices

%%
%set related parameters
alpha = 0.1;
beta = 0.1;
weight = 0.1;
rank = 100;
%run fascinate to get latent feature matrices
[ F ] = fascinate( G,G_new,DO,alpha,beta,weight, rank );
%retore cross-layer dependency matrices with laten feature matrices
D_infer = restoreD(G_new,F);

%generate the connections of newly added nodes
s = double(rand(1,size(F{1}.F,1))>0.8);
% add new node to the corresponding layer
A1 = [G{1}.A;s];
A1 = [A1,[s,0]'];
G{1}.A = A1;
node = size(A1,1);
%run fascinate-zero to find the latent feature of newly added node
[f] = fascinate_zero(A1(node,:),F{1}.F,alpha,beta);