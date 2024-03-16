function norm_vals = arcplot_comparison(all_stacks_data, genes_to_plot, sample_plotting_groups, sample_plotting_colors, gene_grouping, data_type, labels, norm_type, norm_sample_num, save_state)
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
    

    ks_data = cell(1, length(genes_to_plot));
    all_dot_data = {}; all_dot_labels = {};
    f1 = figure; norm_vals = []; all_names = {};
    for f = 1:length(gene_grouping)
        f0 = figure; 
        genes_idx = gene_grouping{f};
        for sp = 1:length(genes_idx)
            ap_ax = subplot(length(genes_idx), 3, 3*(sp - 1) + 1); hold on;
            norm_ap_ax = subplot(length(genes_idx), 3, 3*(sp - 1) + 2); hold on;
            barplot_ax = subplot(length(genes_idx), 3, 3*(sp - 1) + 3); hold on;
            gene_name = '';
            sample_label = 0;         
            
            % plot data
            for gp = 1:length(sample_plotting_groups)
                median_vals = []; std_vals = [];
                for s = 1:length(sample_plotting_groups{gp})
                    num_samples = length(sample_plotting_groups{gp});
                    arc_start = 2*pi*sample_label/num_samples + pi; arc_stop = arc_start - 2*pi/num_samples;
                    sample_label = sample_label + 1;
                    sample_num = sample_plotting_groups{gp}(s);
                    cells_in_stack = cells_data{sample_num};
                                       
                    % only count dots if cell is matched in all rounds
                    [i, ~, ~] = find(cells_in_stack == 0);
                    cell_idx = find(~ismember(1:size(cells_in_stack, 1), unique(i)));
                    dots_in_stack = dots_data{sample_num}; dot_areas_in_stack = dot_area_data{sample_num};

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
                    end
                    dots_to_add = dots_to_add(dots_to_add < prctile(dots_to_add, 99.99));
                    axes(ap_ax); hold on;
                    plot_arc(arc_start, arc_stop, 0, 0, log2(median(dots_to_add)), sample_plotting_colors(sample_label, :), [0 0 0]);
   
                    
                    if strcmp(norm_type, 'median')
                        median_vals(end + 1) = median(dots_to_add);
                    elseif strcmp(norm_type, 'mean')
                        median_vals(end + 1) = mean(dots_to_add);
                    end
                    std_vals(end + 1) = std(dots_to_add);
                    
                end
                
                median_vals = median_vals + 1;
                norm_vals(end + 1) = median_vals(norm_sample_num);
                all_names{end + 1} = gene_name;
                 
                for s = 1:length(sample_plotting_groups{gp})
                    % normalized arc plot
                    arc_start = 2*pi*(s - 1)/num_samples + pi; arc_stop = arc_start - 2*pi/num_samples;
                    axes(norm_ap_ax); plot_arc(arc_start, arc_stop, 0, 0, (median_vals(s))/(norm_vals(end)), sample_plotting_colors(s, :), [0 0 0]);
                    axes(barplot_ax); bar((median_vals)/(norm_vals(end)));
                end
                              
                norm_ap_ax.Title.String = gene_name;
                barplot_ax.Title.String = gene_name;
                barplot_ax.XTick = 1:4;
                barplot_ax.XTickLabels = labels;
                
                bound = ceil(max(median_vals)/(norm_vals(end))) + 0.5;
                
                norm_ap_ax.XLim = bound*[-1 1]; 
                norm_ap_ax.YLim = bound*[-1 1]; 
            end
            
            ap_ax.Title.String = gene_name;
            ap_ax.XLim = [-10 10];
            ap_ax.YLim = [-10 10];


        end   
                
        if save_state
            saveas(f0, ['Samples-' + strjoin(string(labels(sample_plotting_groups{:})), '_') + '-Genes-' + strjoin(string(genes_idx), '_') + norm_type + '_sample' + num2str(norm_sample_num) + '_arcplot.svg'])
        end
    end    
    
    figure(f1);
    bar(norm_vals);
    xticks(1:length(norm_vals));
    xticklabels(all_names);
    xtickangle(45);
    title('Normalizing Values');
    if save_state
        saveas(f1, ['Normalizing_vals_', norm_type, '_sample', num2str(norm_sample_num),'_arcplot.svg']);
        save(['Normalizing_vals_', norm_type, '_sample', num2str(norm_sample_num), '.mat'], 'norm_vals');
    end
end