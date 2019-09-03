function [eigvector, eigvalue, Y] = LPP(X, W, options)
% LPP: Locality Preserving Projections
%
%       [eigvector, eigvalue] = LPP(X, W, options)
% 
%             Input:
%               X       - Data matrix. Each row vector of fea is a data point.
%               W       - Affinity matrix. You can either call "constructW"
%                         to construct the W, or construct it by yourself.
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%                            ReducedDim   -  The dimensionality of the
%                                            reduced subspace. If 0,
%                                            all the dimensions will be
%                                            kept. Default is 0.
%                            PCARatio     -  The percentage of principal
%                                            component kept in the PCA
%                                            step. The percentage is
%                                            calculated based on the
%                                            eigenvalue. Default is 1
%                                            (100%, all the non-zero
%                                            eigenvalues will be kept.
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%               eigvalue  - The eigvalue of LPP eigen-problem. sorted from
%                           smallest to largest. 
% 
% 
%       [eigvector, eigvalue, Y] = LPP(X, W, options) 		
%               
%               Y:  The embedding results, Each row vector is a data point.
%                   Y = X*eigvector
%
%
%    Examples:
%
%       fea = rand(50,70);
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 1;
%       W = constructW(fea,options);
%       options.PCARatio = 0.99
%       [eigvector, eigvalue, Y] = LPP(fea, W, options);
%       
%       
%       fea = rand(50,70);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.bLDA = 1;
%       W = constructW(fea,options);      
%       options.PCARatio = 1;
%       [eigvector, eigvalue, Y] = LPP(fea, W, options);
% 
% 
% Note: After applying some simple algebra, the smallest eigenvalue problem:
%				X^T*L*X = \lemda X^T*D*X
%      is equivalent to the largest eigenvalue problem:
%				X^T*W*X = \beta X^T*D*X
%		where L=D-W;  \lemda= 1 - \beta.
% Thus, the smallest eigenvalue problem can be transformed to a largest 
% eigenvalue problem. Such tricks are adopted in this code for the 


if ~isfield(options,'PCARatio')
    [eigvector_PCA, eigvalue_PCA, meanData, new_X] = PCA(X);
else
    PCAoptions = [];
    PCAoptions.PCARatio = options.PCARatio;
    [eigvector_PCA, eigvalue_PCA, meanData, new_X] = PCA(X,PCAoptions);
end
    
old_X = X;
X = new_X;


[nSmp,nFea] = size(X);

if nFea > nSmp%(LPP必须要求特征数大于样本数
    error('X is not of full rank in column!!');
end

if ~isfield(options,'ReducedDim')
    ReducedDim = nFea; 
else
    ReducedDim = options.ReducedDim; 
end

if ReducedDim > nFea
    ReducedDim = nFea; 
end


D = diag(sum(W));
% L = D - W;
L = D*W;

DPrime = X'*D*X;
DPrime = (DPrime+DPrime')/2;
LPrime = X'*L*X;
LPrime = (LPrime+LPrime')/2;    

dimMatrix = size(DPrime,2);

if dimMatrix > 200 & ReducedDim < dimMatrix/2  % using eigs to speed up!
    option = struct('disp',0);
    [eigvector, eigvalue] = eigs(LPrime,DPrime,ReducedDim,'la',option);
    eigvalue = diag(eigvalue);
else
    [eigvector, eigvalue] = eig(LPrime,DPrime);
    eigvalue = diag(eigvalue);
    
    [junk, index] = sort(-eigvalue);
    eigvalue = eigvalue(index);
    eigvector = eigvector(:,index);
end

eigvalue = ones(length(eigvalue),1) - eigvalue;

if ReducedDim < size(eigvector,2)
    eigvector = eigvector(:, 1:ReducedDim);
    eigvalue = eigvalue(1:ReducedDim);
end

for i = 1:size(eigvector,2)
    eigvector(:,i) = eigvector(:,i)./norm(eigvector(:,i));
end

eigvector = eigvector_PCA*eigvector;


if nargout == 3
    Y = old_X * eigvector;
end

