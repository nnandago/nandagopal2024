function all_data = compare_samples(root_filename, sample_data_filenames, frames_to_analyze)
    
    %root_filename = pwd;
    %sample_data_filenames = {'noBMP_highthres_w_intensities.mat', 'BMP_1h_highthres_w_intensities.mat'};
    num_samples = length(sample_data_filenames);
    stacks = cell(1, num_samples);
    all_data = struct;
    
    all_data.cell_frames_data = cell(1, num_samples);
    all_data.cell_areas_data = cell(1, num_samples);
    all_data.cells_data = cell(1, num_samples);
    
    all_data.dots_data = cell(1, num_samples);
    all_data.dot_mean_intensity_data = cell(1, num_samples);
    all_data.dot_max_intensity_data = cell(1, num_samples);
    all_data.dot_area_data = cell(1, num_samples);
    
    all_data.bg_dots_data = cell(1, num_samples);
    all_data.bg_dot_mean_intensity_data = cell(1, num_samples);
    all_data.bg_dot_max_intensity_data = cell(1, num_samples);
    all_data.bg_dot_area_data = cell(1, num_samples);
    
    all_data.gene_labels = cell(1, num_samples);
    for k = 1:num_samples
        if isfile(fullfile(root_filename, sample_data_filenames{k}))
           display(sample_data_filenames{k})
           var = load(fullfile(root_filename, sample_data_filenames{k}));
           sample_stack = var.stack;
           if isfield(sample_stack.seg_parameters, 'window')
               delete(sample_stack.seg_parameters.window);
           end
           stacks{k} = sample_stack;
           [cell_frame_ids, cells_in_stack, cell_areas_in_stack, dots_in_stack, dot_mean_intensities, dot_max_intensities, dot_areas, bg_dots_in_stack, bg_dot_mean_intensities, bg_dot_max_intensities, bg_dot_areas, gene_labels] ...
               = tabulate_data(sample_stack, frames_to_analyze{k});
           all_data.cell_frames_data{k} = cell_frame_ids;
           all_data.cell_areas_data{k} = cell_areas_in_stack;
           all_data.cells_data{k} = cells_in_stack;
           all_data.dots_data{k} = dots_in_stack;
           all_data.dot_mean_intensity_data{k} = dot_mean_intensities;
           all_data.dot_max_intensity_data{k} = dot_max_intensities;
           [cum_dot_intensity, max_intensity_per_cell, cum_dot_area] = calculate_cumulative_dot_values(dots_in_stack, dot_mean_intensities, dot_max_intensities, dot_areas);
           all_data.dot_area_data{k} = cum_dot_area;
           all_data.dot_cum_intensity_data{k} = cum_dot_intensity;
           all_data.max_intensity_per_cell{k} = max_intensity_per_cell;
           
           all_data.bg_dots_data{k} = bg_dots_in_stack;
           all_data.bg_dot_mean_intensity_data{k} = bg_dot_mean_intensities;
           all_data.bg_dot_max_intensity_data{k} = bg_dot_max_intensities;
           all_data.bg_dot_area_data{k} = bg_dot_areas;
           
           all_data.gene_labels{k} = gene_labels; 
        else
            error(['Cannot open file ', sample_data_filenames{k}]);
            %return;
        end
    end
end





