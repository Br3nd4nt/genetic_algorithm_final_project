clc;
clear all;

addpath('../P1-P12/')

open("P1.mat");

% test
model = ans.ans;
numberOfRuns = 2;
maxIt = 1; 
nPop = 100;


% results = [];
GA_res = GA(model, numberOfRuns, maxIt, nPop);
GEA_1_res = GEA_1(model, numberOfRuns, maxIt, nPop);
GEA_2_res = GEA_2(model, numberOfRuns, maxIt, nPop);
GEA_3_res = GEA_3(model, numberOfRuns, maxIt, nPop);
GEA_res = GEA(model, numberOfRuns, maxIt, nPop);
GEA_tabu_res = GEA_tabu(model, numberOfRuns, maxIt, nPop);
results = [GA_res GEA_1_res GEA_2_res GEA_3_res GEA_res GEA_tabu_res];
