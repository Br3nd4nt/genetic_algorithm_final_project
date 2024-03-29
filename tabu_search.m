classdef tabu_search
methods
    function neighbours = get_neighbours(obj, solution)
        epsilon = 0.05;
        
        new_sol.Position=solution.Position(:, :);
        new_sol.Cost=[];
        new_sol.Rank=[];
        new_sol.DominationSet=[];
        new_sol.DominatedCount=[];
        new_sol.CrowdingDistance=[];
        new_sol.age = 0;
        
        new_sol.age = 0;
        new_sol.operations = 0;
        
        neighbours = repmat(new_sol, size(solution.Position, 2) * 2, 1);
        
        for i = 1:size(solution.Position, 2)
            neighbours(2 * i - 1, 1).Position(1, i) = max(neighbours(2 * i - 1, 1).Position(1, i) - epsilon, 0);
            neighbours(2 * i, 1).Position(1, i) = min(neighbours(2 * i, 1).Position(1, i) + epsilon, 1);
        end
    end

    function best_solution = search(obj, model, init_solution, max_iter, tabu_list_size) 
        % model=model();
        
        CostFunction=@(x) MyCost(x,model);
        
        best_solution = init_solution;
        current_solution = init_solution;
        

        % fix later
        new_sol.Position=[];
        new_sol.Cost=[];
        new_sol.Rank=[];
        new_sol.DominationSet=[];
        new_sol.DominatedCount=[];
        new_sol.CrowdingDistance=[];
        new_sol.age = 0;
        new_sol.operations = 0;
        tabu_list = repmat(new_sol, 1, 1);

        for ignored=1:max_iter
            neighbours = get_neighbours(obj, current_solution);
            best_neighbour = new_sol;
            best_neighbour_fitness = 1e20;
            
            for k = 1:size(neighbours, 1)
                % disp(init_solution.Position)
                neighbour = neighbours(k, 1);
                % if ismember(neighbour, tabu_list) == 0
                if ~any(arrayfun(@(x) isequal(neighbour.Position, x.Position), tabu_list)) 
                    if neighbour.Position ~= init_solution.Position
                        neighbour_fitness =  CostFunction(neighbour.Position);
                        neighbour.Cost = neighbour_fitness;
                        if neighbour_fitness(1) < best_neighbour_fitness
                            best_neighbour_fitness = neighbour_fitness;
                            best_neighbour = neighbour;
                        end
                    end
                end
            end

            % if best_neighbour == new_sol
            if isequaln(best_neighbour, new_sol)
               return;
            end

            current_solution = best_neighbour;
            % disp(best_neighbour)
            % disp(tabu_list(1, 1))
            tabu_list = [tabu_list
                best_neighbour];

            if size(tabu_list, 1) > tabu_list_size
                tabu_list(1) = [];
            end
            if CostFunction(best_neighbour.Position) > CostFunction(best_solution.Position) 
                best_solution = best_neighbour;
            
            end
        end
    end
end
end