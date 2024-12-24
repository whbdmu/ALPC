clear
clc
addpath('data\');
addpath('functions\');
addpath('measure\');
%% dataset 0.6043, 0.5531, 0.6413, 0.5155
load("Wiki_fea.mat");
k = length(unique(Y));
n = length(Y);
for i = 1:length(X)
    X{i} = X{i}';
end
%% param setting
m = 8;
alpha = 1e4;%1 2 3 4
beta = 1e-1;
max = 100;
[A, Z,obj,iter] = ALPC(X, max, k, m, alpha,beta);
[res] = myNMIACCwithmean(Z',Y,k); 
fprintf("ACC, NMI, Purity, F: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f ",  res(1), res(2),res(3),res(4));
