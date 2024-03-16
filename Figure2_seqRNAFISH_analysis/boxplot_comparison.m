function [ks_data, all_dot_data, all_dot_labels]  = boxplot_comparison(all_stacks_data, genes_to_plot, sample_plotting_groups, sample_plotting_colors, gene_grouping, data_type, labels, save_state)
%%% frames_to_plot is a cell of lists, one per filename
%%% genes_to_plot is a cell of gene names or [round, channel] combinations
%%% grouping_per_plot is a cell of indices lists, each being used to group
%%% box plots
%%% gene_grouping is a cell of indices lists, each index list corresponding
%%% to genes that should be shown as subplots within a single figure

    % all_stacks_data = compare_samples(root_dir, stack_filenames, frames_to_plot);
    
    % cell array of list of ids within each stack that each row in subsequent data
    % corresponds to
    cells_data = all_stacks_data.cells_data;
    
    % cell array of list of frames within the stack that each row in subsequent data
    % corresponds to
    cell_frames_data = all_stacks_data.cell_frames_data;
    
    % cell array of 2D matrices of cell areas -> rows are matched cells, columns
    % are rounds
    cell_areas_data = all_stacks_data.cell_areas_data;
    
    % cell array of 3D matrices of number of dots -> rows are matched cells,
    % columns are round, depth is channels
    dots_data = all_stacks_data.dots_data;
    
    % cell arrayz of 3D cell matrices with each element containing a list of dot intensities
    dot_mean_intensity_data = all_stacks_data.dot_mean_intensity_data;
    dot_max_intensity_data = all_stacks_data.dot_max_intensity_data;
    dot_area_data = all_stacks_data.dot_area_data;
    dot_cum_intensity_data = all_stacks_data.dot_cum_intensity_data;
    if isfield(all_stacks_data, 'nuc_signal_sum')
        nuc_signal_data =  all_stacks_data.nuc_signal_sum;
    end
    
    num_samples = length(cell_frames_data);

    ks_data = cell(1, length(genes_to_plot));
    all_dot_data = {}; all_dot_labels = {};
    for f = 1:length(gene_grouping)
        f0 = figure;
        genes_idx = gene_grouping{f};
        for sp = 1:length(genes_idx)
            bp_ax = subplot(length(genes_idx), 3, 3*(sp - 1) + 1);
            consolidated_dot_data = [];
            grouping_vars_level1 = [];
            grouping_vars_level2 = [];
            gene_name = '';
            sample_label = 0;         
            
            % plot data
            for gp = 1:length(sample_plotting_groups)
                for s = 1:length(sample_plotting_groups{gp})
                    sample_label = sample_label + 1;
                    sample_num = sample_plotting_groups{gp}(s);
                    cells_in_stack = cells_data{sample_num};
                    
                    
                    % only count dots if cell is matched in all rounds
                    [i, ~, ~] = find(cells_in_stack == 0);
                    cell_idx = find(~ismember(1:size(cells_in_stack, 1), unique(i)));
                    dots_in_stack = dots_data{sample_num}; dot_areas_in_stack = dot_area_data{sample_num};
                    
                    % filter out sparse cells
%                     cum_dots_in_round = [];
%                     for r = 1:5
%                         cum_dots_in_round = cat(2, cum_dots_in_round, sum(squeeze(dots_in_stack(:, r, 3:5)), 2));
%                     end
%                     
%                    
%                     avg_dots_per_cell = sum(cum_dots_in_round, 2); %sum(sum(dots_in_stack, 3), 2)/14;
%                     sparse_cells_idx = find(avg_dots_per_cell < 200);
%                     cell_idx = cell_idx(~ismember(cell_idx, sparse_cells_idx));
                    

                    % find [r, c] index of the gene within the stack, if
                    % necessary
                    if ischar(genes_to_plot{1})
                        gene_labels = all_stacks_data.gene_labels{1};
                        gene_name = genes_to_plot{genes_idx(sp)};
                        if ~isempty(strcmp(gene_name, gene_labels))
                            [r, c] = ind2sub([size(dots_in_stack, 2), size(dots_in_stack, 3)], find(strcmp(gene_name, gene_labels)));
                        else
                            return;
                        end
                    else
                        gene_name = [num2str(genes_to_plot{genes_idx(sp)}(1)), '-', num2str(genes_to_plot{genes_idx(sp)}(2))];
                    end

                    if strcmp(data_type, 'count')
                        dots_to_add = dots_in_stack(cell_idx, r, c);
                    elseif strcmp(data_type, 'area')
                        dots_to_add = dot_areas_in_stack(cell_idx, r, c);
                    elseif strcmp(data_type, 'intensity')
                        dots_to_add = dot_cum_intensity_data{sample_num}(cell_idx, r, c);
                    elseif strcmp(data_type, 'nuc')
                        dots_to_add = nuc_signal_data{sample_num}(cell_idx, r, c);
                    elseif strcmp(data_type, 'cov')
                        dots_to_add = dots_in_stack(cell_idx, r, c);
                    end
                    dots_to_add = dots_to_add(dots_to_add < prctile(dots_to_add, 99.9)); % & dots_to_add > prctile(dots_to_add, 0.05));
                    consolidated_dot_data = vertcat(consolidated_dot_data, dots_to_add);
                    grouping_vars_level1 = horzcat(grouping_vars_level1, sample_label*ones(1, length(dots_to_add)));
                end
                grouping_vars_level2 =  horzcat(grouping_vars_level2, gp*ones(1, length(grouping_vars_level1) - length(grouping_vars_level2)));
            end
            
            grouping_labels = vertcat(grouping_vars_level1, grouping_vars_level2)';
            all_dot_data{end + 1} = consolidated_dot_data; all_dot_labels{end + 1} = grouping_vars_level1;
   
            if strcmp(data_type, 'cov')
                cov_vals = zeros(sample_label, 1);
                for s = 1:sample_label
                    cov_vals(s) = std(consolidated_dot_data(grouping_vars_level1 == s))/mean(consolidated_dot_data(grouping_vars_level1 == s));
                end
                bar(1:sample_label, cov_vals);
                ylabel('Coeff. of Variation');
                xticklabels(labels(sample_plotting_groups{:}));
                ylim([0 1.5]);
                title(gene_name);
            else
                boxplot(bp_ax, consolidated_dot_data, grouping_labels, 'Color', sample_plotting_colors, 'FactorGap', [3 1]);
                bp_ax.Title.String = gene_name;
                bp_ax.YLabel.String = '# Transcripts';
                if prctile(consolidated_dot_data, 1) ~= prctile(consolidated_dot_data, 95)
                    bp_ax.YLim = [prctile(consolidated_dot_data, 1) prctile(consolidated_dot_data, 99)];
                else
                    % do nothing;
                end
                bp_ax.XTickLabel = labels(sample_plotting_groups{:});
            end
            
%             vp_ax = subplot(length(genes_idx), 2, 3*(sp - 1) + 2 );
             vp_ax = subplot(length(genes_idx), 3, 3*(sp - 1) + 2);
            if ~isempty(sample_plotting_colors)
                violinplot(consolidated_dot_data, grouping_vars_level1, 'ShowMean', true, 'ViolinColor', sample_plotting_colors);
            else
                violinplot(consolidated_dot_data, grouping_vars_level1, 'ShowMean', true);
            end
            vp_ax.Title.String = gene_name;
            vp_ax.YLabel.String = '# Transcripts ';
            vp_ax.YLim = [prctile(consolidated_dot_data, 1) prctile(consolidated_dot_data, 99)];
            vp_ax.XTickLabel = labels(sample_plotting_groups{:});
            
            vp_log_ax = subplot(length(genes_idx), 3, 3*sp);
            if ~isempty(sample_plotting_colors)
                violinplot(log2(consolidated_dot_data + 1), grouping_vars_level1, 'ShowMean', true, 'ViolinColor', sample_plotting_colors);
            else
                violinplot(log2(consolidated_dot_data + 1), grouping_vars_level1, 'ShowMean', true);
            end
            vp_log_ax.Title.String = gene_name;
            vp_log_ax.YLabel.String = '# Transcripts (log2)';
            vp_log_ax.YLim = [prctile(log2(consolidated_dot_data + 1), 1) prctile(log2(consolidated_dot_data), 99.99)];
            vp_log_ax.XTickLabel = labels(sample_plotting_groups{:});
            
            
%             cov_ax = subplot(length(genes_idx), 3, 3*sp);
%             cov_vals = zeros(sample_label, 1);
%             for s = 1:sample_label
%                 cov_vals(s) = std(consolidated_dot_data(grouping_vars_level1 == s))/mean(consolidated_dot_data(grouping_vars_level1 == s));
%             end
%             bar(1:sample_label, cov_vals);
            
            % calculate KS significance within plotting groups
            gene_ks_data = struct;
            unique_sample_ids = unique(grouping_vars_level1);
            h_data = zeros(length(unique_sample_ids), length(unique_sample_ids));
            p_data = zeros(length(unique_sample_ids), length(unique_sample_ids));
 
            for s1 = 1:length(unique_sample_ids) - 1
                for s2 = (s1 + 1):length(unique_sample_ids)
                    [h, p] = kstest2(consolidated_dot_data(grouping_vars_level1 == unique_sample_ids(s1)), consolidated_dot_data(grouping_vars_level1 == unique_sample_ids(s2)));
                    h_data(s1, s2) = h; h_data(s2, s1) = h;
                    p_data(s1, s2) = p; p_data(s2, s1) = p;
                end
            end
            gene_ks_data.h_data = h_data;
            gene_ks_data.p_data = p_data;
            ks_data{genes_idx(sp)} = gene_ks_data;
        end   
        if save_state
            saveas(f0, ['Samples-' + strjoin(string(labels(sample_plotting_groups{:})), '_') + '-Genes-' + strjoin(string(genes_idx), '_') + '.svg'])
            if strcmp(data_type, 'cov')
                saveas(f0, ['Samples-' + strjoin(string(labels(sample_plotting_groups{:})), '_') + '-Genes-' + strjoin(string(genes_idx), '_') + 'withCOV.svg'])
            end
        end
        
    end            
end