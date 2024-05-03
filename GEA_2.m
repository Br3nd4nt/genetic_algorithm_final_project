function RESULTS = GEA_2(model, numberOfRuns, maxIt, nPop)
    RESULTS = [];
    runsCount = 0;
    
    CostFunction = @(x) MyCost(x,model);
    nVar = model.d + model.n + model.m + model.p + model.r + model.s + model.j;
    VarSize = [1 nVar];
    VarMin = 0;
    VarMax=1;

    pCrossover=0.8;                             % Crossover Percentage
    nCrossover = 2*round(pCrossover*nPop/2);    % Number of Parnets (Offsprings)
    
    pMutation = 0.2;                            % Mutation Percentage
    nMutation = round(pMutation*nPop);          % Number of Mutants
    
    mu=0.02;                                    % Mutation Rate
    
    sigma = 0.1*(VarMax-VarMin);                % Mutation Step Size
    
    % masks hyperparameters
    masks_bitmap = [1 1 1];
    operations_bitmap = [0 1 0];
    inejction_rate = 0.5;
    epsilon = 0.05;
    mask_p_percent = 0.75;
    single_threshold = 20;
    double_threshold = 14;
    triple_threshold = 10;
    masks = masks_class();

    % mortality 
    max_age = 5;
    max_operations = 5;
    
    % Initialization

    empty_individual.Position=[];
    empty_individual.Cost=[];
    empty_individual.Rank=[];
    empty_individual.DominationSet=[];
    empty_individual.DominatedCount=[];
    empty_individual.CrowdingDistance=[];
    
    empty_individual.age = 0;
    empty_individual.operations = 0;

    while runsCount < numberOfRuns
    runsCount = runsCount + 1;

    % setup
    result.GA       = [];
    result.GEA_1    = [];
    result.GEA_2    = [];
    result.GEA_3    = [];
    result.GEA      = [];
    result.GEA_tabu = [];
    result.result = 0;
    result.cpuTime = 0;
    
    % starting the time
    tic;

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
    
    
    % MAIN LOOP
    for it = 1:maxIt
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
        
        % Mask creation
        masks_to_choose = zeros(nPop * mask_p_percent, size(pop(1).Position, 2), 3);
        if masks_bitmap(1) == 1
        masks_to_choose(:, :, 1) = masks.create_single_gene_mask(pop, epsilon, mask_p_percent, single_threshold);
    
        end
        if masks_bitmap(2) == 1
        masks_to_choose(:, :, 2) = masks.create_double_gene_mask(pop, epsilon, mask_p_percent, double_threshold);
        end
        if masks_bitmap(3) == 1
        masks_to_choose(:, :, 3) = masks.create_triple_gene_mask(pop, epsilon, mask_p_percent, triple_threshold);
        end
        
        % Second operation
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
    
        % Merge
        pop=[pop
             popc
             popm
             pop_sc_2
             ];
    
        % Clearing empty solutions and check ages
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
        

        F1=pop(F{1});
        result.GEA_2=[result.GEA_2; F1(1).Cost(1)];
        disp(['(GEA_2, run: ' num2str(runsCount) ') Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1)) '; Best cost ' num2str(result.GEA_2(it))]);
    end

    result.cpuTime = toc;
    result.result = F1(1).Cost(1);



    % RESULTS(runsCount) = result;
    RESULTS = [RESULTS; result];
    end
end