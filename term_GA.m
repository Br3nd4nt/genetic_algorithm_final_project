% clc;
% clear all;
% close all;
% model=model();
A_CODE_RESULTS = [];
A_CURRENT_ITERATIONS_COUNT = 0;
A_MAX_CODE_ITERATIONS_ = 5;
MaxIt= 200;      % Maximum Number of Iterations
 
while A_CURRENT_ITERATIONS_COUNT < A_MAX_CODE_ITERATIONS_
    A_CURRENT_ITERATIONS_COUNT = A_CURRENT_ITERATIONS_COUNT + 1;

tic
%% Problem Definition

% model=model();

CostFunction=@(x) MyCost(x,model);      % Cost Function

nVar=model.d+model.n+model.m+model.p+model.r+model.s+model.j;   % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables

% Number of Objective Functions
nObj=2;   %Both of them are minimization 


%% NSGA-II Parameters

% MaxIt=1000;      % Maximum Number of Iterations

nPop=100;        % Population Size

pCrossover=0.7;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation=0.2;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

mu=0.02;                    % Mutation Rate

sigma=0.1*(VarMax-VarMin);  % Mutation Step Size


%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
end

GA=zeros(MaxIt,1);

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);


%% NSGA-II Main Loop

for it=1:MaxIt
    
    % Crossover
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        
        i1=randi([1 nPop]);
        p1=pop(i1);
        
        i2=randi([1 nPop]);
        p2=pop(i2);
        
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        
        i=randi([1 nPop]);
        p=pop(i);
        
        popm(k).Position=Mutate(p.Position,mu,sigma);
        
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
    % Merge
    pop=[pop
         popc
         popm]; %#ok
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
    
    % Store F1
    F1=pop(F{1});
    
    GA(it)=F1(1).Cost(1);

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1)) ' Best cost ' num2str(GA(it))]);
    
    % Plot F1 Costs
    % figure(1);
    % PlotCosts(F1);
    % pause(0.01);
    
end
%% metrics
cputime=toc;
mid=MID(F1);
nps=length(F1);
sns=SNS(mid,F1);
ms=MS(F1);
res=mid/ms;
 disp([ ' NPS= '  num2str(nps) '   CPU TIME= ' num2str( cputime)  '   MID= ' num2str(mid) '   SNS=' num2str(sns) '  MS=' num2str(ms) '  Response=' num2str(res) ]); 

%% Results
disp('====================');
z1=pop(1).Cost(1);
z1
% z2=pop(1).Cost(2);
% z2


A_CODE_RESULTS = [A_CODE_RESULTS; pop(1).Cost(1) toc];
end
% plot(GA)