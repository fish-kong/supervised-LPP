function W = constructW3(fea,gnd);

options.k = 20;

options.t = 1;

options.bSelfConnected = 1;

delta=100;

[nSmp, nFea] = size(fea);

D = zeros(nSmp);

%�ල��ŷʽ����
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

% �Ǽලŷʽ����
% for i=1:nSmp-1
%     for j=i+1:nSmp
%      D(i,j)=norm(fea(i,:)-fea(j,:));
%     end
% end


D = D+D';

%������������ͼ
G = zeros(nSmp,nSmp);
[~, idx] = sort(D, 2); % sort each row ����Ϊ��λ ��������
for i=1:nSmp
    G(i,idx(i,1:options.k+1)) = 1;%������i�����k+1������Ϊ1
    % ��������֮��Ϊ1 �ǽ���Ϊ0
end


if ~options.bSelfConnected
    G  = G - diag(diag(G));
end
%������������Ȩ�ؾ���W
D2 = exp(-D.^2/options.t);%����ÿ�����������ƶ�
W = D2.*G;% ��������֮��Ϊ���ƶ�ΪD �ǽ���Ϊ0
% W = max(W,W');
% W = sparse(W);



