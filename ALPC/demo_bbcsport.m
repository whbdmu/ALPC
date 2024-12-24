clear
clc
addpath('data\');
addpath('functions\');
addpath('measure\');
%% dataset 0.9567, 0.8474, 0.9567, 0.9033
load("BBCSport.mat");
k = length(unique(Y));
n = length(Y);
%% param setting
m = 15;
alpha = 90;%1 2 3 4
beta = 1e-1;
max = 100;
[A, Z,obj,iter] = ALPC(X, max, k, m, alpha,beta);
[res] = myNMIACCwithmean(Z',Y,k); 
fprintf("ACC, NMI, Purity, F: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f ",  res(1), res(2),res(3),res(4));
