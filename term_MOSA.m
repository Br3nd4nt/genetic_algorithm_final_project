%% SA 
clc;
% clear all;
close all;
%% Porblem Definition
tic;
model=model(); 

CostFunction=@(x) MyCost(x,model);      % Cost Function

nVar=model.d+model.n+model.m+model.p+model.r+model.s+model.j;   % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables

% Number of Objective Functions
nObj=2;     % Both of them are minimization 
%% SA Parameters 

Maxit=1000;      % Maximum Number of Iterations
T=500;          % Initial temprature 
damp=0.99;     % Rate of reduction 

mu=0.7;                    % Mutation Rate

sigma=0.1*(VarMax-VarMin);  % Mutation Step Size

% Empty 
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];

Sol=repmat(empty_individual,1);
NewSol=repmat(empty_individual,1);
% Initial solution 
Sol.Position=unifrnd(VarMin,VarMax,VarSize);
    
Sol.Cost=CostFunction(Sol.Position);

% ParetoCost=repmat(empty_individual,Maxit,1); 
BestSol=Sol;
SA=zeros(Maxit,1);
% BestCost2=zeros(Maxit,1);   
%% main loop
for it=1:Maxit
       
        NewSol.Position=Mutate(Sol.Position,mu,sigma);
        
        NewSol.Cost=CostFunction(NewSol.Position);
    
     
     if (NewSol.Cost(1)<=Sol.Cost(1))% && (NewSol.Cost(2)<=Sol.Cost(2))
       Sol=NewSol;
       % ParetoCost(it)=Sol; 
     % elseif (NewSol.Cost(1)<=Sol.Cost(1))% || (NewSol.Cost(2)<=Sol.Cost(2))
         % ParetoCost(it)=NewSol;
     elseif (NewSol.Cost(1)>Sol.Cost(1))% && (NewSol.Cost(2)>Sol.Cost(2))         
           delta1=NewSol.Cost(1)-Sol.Cost(1);           
           p1=exp(-delta1/T);
           h=rand;    
           %delta2=NewSol.Cost(2)-Sol.Cost(2);
            %p2=exp(-delta2/T);
           if(p1<=h)%&&(p2<=h)
                  Sol=NewSol;                   
           end
           %ParetoCost(it)=Sol; 
           
     end
     
      if(Sol.Cost(1) <=BestSol.Cost(1))%&&(Sol.Cost(2) <=BestSol.Cost(2))
           BestSol=Sol;
       end  
    SA(it)=BestSol.Cost(1); 
    %BestCost2(it)=BestSol.Cost(2);   
    T=damp*T;
    

    disp(['Iteration' num2str(it) ' : Best Cost1 =' num2str(SA(it))]);  % ' : Best Cost2 =' num2str(BestCost2(it))]);  
    
end  
%% Pareto fronts 
% disp('=Please waite; Comouter is computing Pareto fronts=');
%     % Non-Dominated Sorting
%     [pop, F]=NonDominatedSorting(ParetoCost);
% 
%     % Calculate Crowding Distance
%     pop=CalcCrowdingDistance(pop,F);
% 
%     % Sort Population
%     pop=SortPopulation(pop);
% 
%     % Non-Dominated Sorting
%     [pop, F]=NonDominatedSorting(pop);
% 
%     % Calculate Crowding Distance
%     pop=CalcCrowdingDistance(pop,F);
% 
%     % Sort Population
%     [pop, F]=SortPopulation(pop);
% 
%     % Store F1
%     F1=pop(F{1});
% 
%     % Show Iteration Information
%     disp([ ': Number of F1 Members = ' num2str(numel(F1))]);
    
%% metric
% cputime=toc;
% mid=MID(F1);
% nps=length(F1);
% sns=SNS(mid,F1);
% ms=MS(F1);
% res=mid/ms;
%  disp([ ' NPS= '  num2str(nps) '   CPU TIME= ' num2str( cputime)  '   MID= ' num2str(mid) '   SNS=' num2str(sns) '  MS=' num2str(ms) '  Response=' num2str(res) ]); 


 %% Results
disp('====================');
disp(BestSol)
disp('====================');
disp("CPU TIME")
disp(toc)
% z1=BestSol.Cost%(1);
% z1
% z2=BestSol.Cost(2);
% z2
%% figures
% figure(1);
% plot(SA);
% figure(2);
% plot(BestCost2);
% figure(3);
% PlotCosts(F1);

