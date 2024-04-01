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


%% GEA Parameters

% MaxIt=1000;      % Maximum Number of Iterations

nPop=100;        % Population Size

pCrossover=0.7;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation=0.2;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

mu=0.02;                    % Mutation Rate

sigma=0.1*(VarMax-VarMin);  % Mutation Step Size

%% masks hyperparameters
masks_bitmap = [1 1 1];
operations_bitmap = [0 1 0];
GEA_2=zeros(MaxIt,1);
inejction_rate = 0.4;
epsilon = 0.05;
mask_p_percent = 0.6;
single_threshold = 12;
double_threshold = 10;
triple_threshold = 8;
masks = masks_class();

%% mortality 
max_age = 20;
max_operations = 10;

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];

%% FOR GEA ONLY
empty_individual.age = 0;
empty_individual.operations = 0;


pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
end



% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);


%% NSGA-II Main Loop

for it=1:MaxIt
    
    for len=1:size(pop, 1) 
        pop(len).age = pop(len).age + 1;
    end

    % Crossover
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        
        i1=randi([1 nPop]);
        p1=pop(i1);

        i2=randi([1 nPop]);
        p2=pop(i2);

        pop(i1).operations = pop(i1).operations + 1;
        pop(i2).operations = pop(i2).operations + 1;
        
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
        pop(i).operations = pop(i).operations + 1;
        
        popm(k).Position=Mutate(p.Position,mu,sigma);
        
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
    %mask creation
    masks_to_choose = zeros(100 * mask_p_percent, size(pop(1).Position, 2), 3);
    if masks_bitmap(1) == 1
    masks_to_choose(:, :, 1) = masks.create_single_gene_mask(pop, epsilon, mask_p_percent, single_threshold);

    end
    if masks_bitmap(2) == 1
    masks_to_choose(:, :, 2) = masks.create_double_gene_mask(pop, epsilon, mask_p_percent, double_threshold);
    end
    if masks_bitmap(3) == 1
    masks_to_choose(:, :, 3) = masks.create_triple_gene_mask(pop, epsilon, mask_p_percent, triple_threshold);
    end

    % first operation
    if operations_bitmap(1) == 1
    pop_sc_1_vals = masks.first_scenario(pop, masks.get_random_mask(masks_bitmap, masks_to_choose), pCrossover);
    pop_sc_1=repmat(empty_individual,size(pop_sc_1_vals, 1),1);

    for k=1:size(pop_sc_1_vals, 1)
        pop_sc_1(k).Position = pop_sc_1_vals(k, :);

        pop_sc_1(k).Cost=CostFunction(pop_sc_1(k).Position);
    end
    else
        pop_sc_1 = [];
    end
    
    %second operation
    if operations_bitmap(2) == 1
        pop_sc_2 = repmat(empty_individual,nMutation,1);
    

        for k=1:nMutation
            
            i=randi([1 nPop]);
            p=pop(i);
            pop(i).operations = pop(i).operations + 1;
            pop_sc_2(k).Position=masks.second_scenario(p.Position, mu, sigma, masks.get_random_mask(masks_bitmap, masks_to_choose));
            
            pop_sc_2(k).Cost=CostFunction(pop_sc_2(k).Position);
            
        end
    else 
        pop_sc_2 = [];
    end

    %third operaion
    if operations_bitmap(3) == 1
        pop_sc_3 = masks.third_scenario(pop, masks.get_random_mask(masks_bitmap, masks_to_choose), inejction_rate);
    
        for k=1:size(pop_sc_3, 1)
    
            pop_sc_3(k).Cost=CostFunction(pop_sc_3(k).Position);
        end
    else 
        pop_sc_3= [];
    end
    
    % Merge
    pop=[pop
         popc
         popm
         pop_sc_1
         pop_sc_2
         pop_sc_3
         ]; %#ok

% clearing empty solutions and check ages
k = numel(pop);
while k >= 1
    if k > numel(pop) 
        k = k - 1;
        continue
    else 
        if isempty(pop(k).Cost) ||  pop(k).age > max_age || pop(k).operations > max_operations
            pop(k) = [];
        end
    end
    k = k - 1;
end

    for k=1:numel(pop)
        if size(pop(k).Position, 1) > 1
            pop(k).Position = pop(k).Position';
        end
    end



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
    
    GEA_2(it)=F1(1).Cost(1);

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1)) ' Best cost ' num2str(GEA_2(it))]);
    % disp(['Iteration ' num2str(GEA(it))]);
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
% z1
A_CODE_RESULTS = [A_CODE_RESULTS; pop(1).Cost(1) toc];
end
% plot(GEA)