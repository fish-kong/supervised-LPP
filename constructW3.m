function W = constructW3(fea,gnd);

options.k = 20;

options.t = 1;

options.bSelfConnected = 1;

delta=100;

[nSmp, nFea] = size(fea);

D = zeros(nSmp);

%监督核欧式距离
for i=1:nSmp-1
    for j=i+1:nSmp
        d=sqrt(2-2*exp(-norm(fea(i,:)-fea(j,:)).^2/delta));
        
        if gnd(i)==gnd(j)
            D(i,j) = sqrt(1-exp(-d^2/10));
        else
            D(i,j) =sqrt( exp(-d^2/10));
        end
    end
end

% 非监督欧式距离
% for i=1:nSmp-1
%     for j=i+1:nSmp
%      D(i,j)=norm(fea(i,:)-fea(j,:));
%     end
% end


D = D+D';

%构建近邻无向图
G = zeros(nSmp,nSmp);
[~, idx] = sort(D, 2); % sort each row 以行为单位 升序排列
for i=1:nSmp
    G(i,idx(i,1:options.k+1)) = 1;%和样本i最近的k+1个样本为1
    % 近邻样本之间为1 非近邻为0
end


if ~options.bSelfConnected
    G  = G - diag(diag(G));
end
%构建近邻链接权重矩阵W
D2 = exp(-D.^2/options.t);%计算每个样本的相似度
W = D2.*G;% 近邻样本之间为相似度为D 非近邻为0
% W = max(W,W');
% W = sparse(W);



