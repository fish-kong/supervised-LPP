clc;clear;close all
load wine;fea=double(wine);gnd=wine_labels;

%% PCA
% options = [];
% options.ReducedDim = 3;
% options.PCARatio = 1;
% [eigvector,eigvalue,meanData, Y] = PCA(fea,options);

%%
W = constructW3(fea,gnd);

options.ReducedDim = 3;
options.PCARatio = 1;
[eigvector, eigvalue, Y] = LPP(fea, W, options);
%% 

%%
mappedX1=Y;
figure
for i=1:max(gnd)
    n=find(gnd==i);
    plot(mappedX1(n,1),mappedX1(n,2),'*')
    hold on
    grid on
end
figure

for i=1:max(gnd)
    n=find(gnd==i);
    plot3(mappedX1(n,1),mappedX1(n,2),mappedX1(n,3),'*')
    hold on
    grid on
end
legend('1','2','3')
