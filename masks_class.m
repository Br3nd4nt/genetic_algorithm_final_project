classdef masks_class
methods
    function mask = get_random_mask(~, masks_bitmap, masks)
            randomIndex = randi(numel(find(masks_bitmap == 1)));
            mask = masks(:, :, randomIndex);
    end

    function mask = create_single_gene_mask(~, population, epsilon, p, threshold)
            rows = size(population, 1);
            cols = size(population(1).Position, 2);
            mask_rows = ceil(rows * p);
            mask = zeros(mask_rows, cols);
            
            for j=1:mask_rows
                for i = 1:cols
                    if (mask(j, i) == 1)
                        continue
                    end
                    counter = 0;
                    indicies = [];
                    for col_iter=1:mask_rows
                        if (abs(population(j).Position(i) - population(col_iter).Position(i)) <= epsilon) 
                            counter = counter + 1;
                            indicies(size(indicies) + 1) = col_iter;
                        end
                    end
                    if counter >= threshold
                        for t=indicies
                            mask(t, i) = 1;
                        end
                    end
                end
            end
    end

    function mask = create_double_gene_mask(~, population, epsilon, p, threshold)
            rows = size(population, 1);
            cols = size(population(1).Position, 2);
            mask_rows = ceil(rows * p);
            mask = zeros(mask_rows, cols);
            
            for j=1:mask_rows
                for i = 1:cols - 1
                    if (mask(j, i) == 1 || mask(j, i + 1) == 1)
                        continue
                    end
                    counter = 0;
                    indicies = [];
                    for col_iter=1:mask_rows
                        if (abs(population(j).Position(i) - population(col_iter).Position(i)) <= epsilon) && (abs(population(j).Position(i + 1) - population(col_iter).Position(i + 1)) <= epsilon) 
                            counter = counter + 1;
                            indicies(size(indicies) + 1) = col_iter;
                        end
                    end
                    if counter >= threshold
                        for t=indicies
                            mask(t, i) = 1;
                            mask(t, i + 1) = 1;
                        end
                    end
                end
            end
    end

    function mask = create_triple_gene_mask(~, population, epsilon, p, threshold)
            rows = size(population, 1);
            cols = size(population(1).Position, 2);
            mask_rows = ceil(rows * p);
            mask = zeros(mask_rows, cols);
            
            for j=1:mask_rows
                for i = 1:(cols - 2)
                    if (mask(j, i) == 1 || mask(j, i + 1) == 1 || mask(j, i + 2) == 1)
                        continue
                    end
                    counter = 0;
                    indicies = [];
                    for col_iter=1:mask_rows
                        if ((abs(population(j).Position(i) - population(col_iter).Position(i)) <= epsilon)) && (abs(population(j).Position(i + 1) - population(col_iter).Position(i + 1)) <= epsilon) && (abs(population(j).Position(i + 2) - population(col_iter).Position(i + 2)) <= epsilon) 
                            counter = counter + 1;
                            indicies(size(indicies) + 1) = col_iter;
                        end
                    end
                    if counter >= threshold
                        for t=indicies
                            mask(t, i) = 1;
                            mask(t, i + 1) = 1;
                            mask(t, i + 2) = 1;
                        end
                    end
                end
            end
    end

    function offsprings = first_scenario(~, population, mask, crossover_rate)
           % find dominant
           vals = sum(mask, 2);
           [~, I] = max(vals);
           population(I).operations = population(I).operations + 1;

           number_of_crossover = ceil(size(population, 1) * crossover_rate / 2);
           offsprings = zeros(number_of_crossover * 2, size(population(1).Position, 2));
           for i = 1:number_of_crossover
                second_parent = randi(size(population, 1));
                while (I == second_parent)
                    second_parent = randi(size(population, 1));
                end
                
                population(second_parent).operations = population(second_parent).operations + 1;

                crossover_point = randi([1, size(population(1).Position, 2)]);

                offsprings(i * 2, 1:crossover_point) = population(I).Position(1, 1:crossover_point);
                offsprings(i * 2, crossover_point:end) = population(second_parent).Position(1, crossover_point:end);
                offsprings(i * 2 + 1, 1:crossover_point) = population(second_parent).Position(1, 1:crossover_point);
                offsprings(i * 2 + 1, crossover_point:end) = population(I).Position(1, crossover_point:end);
           end

    end

    function y = second_scenario(~, x, mu, sigma, mask) % idunno shuold work
        nVar = numel(x);
        
        nMu = ceil(mu * nVar);
        j = randsample(nVar, nMu);
        if numel(sigma) > 1
            sigma = sigma(j);
        end
        % y = x.copy();
        y = x(:);
        for i = 1:nMu
            % disp(mask)
            % disp(nMu)
            if mask(j) == 0  % Only mutate if mask bit is 0
                y(j(i)) = x(j(i)) + (sigma .* x(j(i)));  % Use sigma(i) for individual mutations
            end
        end
    end

    function injected = third_scenario(~, population, mask, injection_rate)
            number_of_injections = ceil((size(population, 1) - size(mask, 1)) * injection_rate);
            % find dominant
            vals = sum(mask, 2);
            [~, I] = max(vals);
            dominant = population(I);
            dominant.operations = dominant.operations + 1;
            dominant_mask = mask(I, :);
            indicies = randperm((size(population, 1) - size(mask, 1)), number_of_injections) + size(mask, 1);
            
            injected = population(indicies, :);
            
            for i=1:size(injected, 1)
                for j=1:size(mask, 2)
                    if dominant_mask(1, j) == 1
                        injected(i).Position(j) = dominant.Position(j);
                    end
                end
            end
        end
end
end