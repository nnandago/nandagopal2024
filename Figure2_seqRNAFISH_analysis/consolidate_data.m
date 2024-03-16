function dot_data_tables = consolidate_data(all_stacks_data, gene_list, data_type, normalize_type)
%%% outputs a cell array of tables, each corresponding to the counts
%%% matrix for a stack

%%% gene list is either a cell array of strings or [r, c] combinations

    cells_data = all_stacks_data.cells_data;
    dots_data = all_stacks_data.dots_data;
    dot_area_data = all_stacks_data.dot_area_data;

    num_samples = length(all_stacks_data.cell_frames_data);
    dot_data_tables = cell(1, num_samples);
    for s = 1:num_samples
        cells_in_stack = cells_data{s};
        
        % only count dots if cell is matched in all rounds
        [i, ~, ~] = find(cells_in_stack == 0);
        cell_idx = find(~ismember(1:size(cells_in_stack, 1), unique(i)));
        dots_in_stack = dots_data{s}; dot_areas_in_stack = dot_area_data{s};

        
        new_table = zeros(length(cell_idx), length(gene_list));
        for gn = 1:length(gene_list)
            if ischar(gene_list{1})
                gene_labels = all_stacks_data.gene_labels{s};
                gene_name = gene_list{gn};
                if ~isempty(strcmp(gene_name, gene_labels))
                    [r, c] = ind2sub([size(dots_in_stack, 2), size(dots_in_stack, 3)], find(strcmp(gene_name, gene_labels)));
                else
                    return;
                end
            else
                r = gene_list{s}(1);
                c = gene_list{s}(2);

            end
            if isempty(dots_in_stack(cell_idx, r, c))
                keyboard;
            end
            
            if strcmp(data_type, 'count')
                dots_to_add = dots_in_stack(cell_idx, r, c);
            elseif strcmp(data_type, 'area')
                dots_to_add = dot_areas_in_stack(cell_idx, r, c);
            end
         
            new_table(:, gn) = dots_to_add;
        end
            
        if s == 1
            if strcmp(normalize_type, 'median')
                norm_vector = median(new_table, 1);
            elseif strcmp(normalize_type, 'mean')
                norm_vector = mean(new_table, 1);
            else
                norm_vector = ones(1, size(new_table, 2));
            end
        end
           
        %throw out outliers
        for col = 1:size(new_table, 2)
            new_table(new_table(:, col) > prctile(new_table(:, col), 99), col) = prctile(new_table(:, col), 99);
        end
        dot_data_tables{s} = array2table(bsxfun(@rdivide, new_table, norm_vector), 'VariableNames', gene_list);
    
    end

end