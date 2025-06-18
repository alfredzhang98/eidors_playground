function eidors_interactive_ui()
% EIDORS Interactive UI - Interactive EIDORS simulation interface
% Fully implemented according to EIDORS official tutorial standards

    mytoolbox_root = fileparts(mfilename('fullpath'));  
    eidors_folder    = fullfile(mytoolbox_root, 'eidors');
    startup_script   = fullfile(eidors_folder, 'eidors_startup.m');
    if exist(startup_script,'file')
    run(startup_script);
    else
    warning('eidors_startup.m cannot find please intall it.');
    end

    % Check EIDORS
    if ~exist('eidors_obj', 'file')
        error('Please initialize EIDORS first! Run startup_eidors command.');
    end
    
    fprintf('Starting EIDORS interactive interface...\n');
    
    % Use persistent variables to avoid parameter passing issues
    persistent app_data;
    
    if isempty(app_data)
        app_data = init_app_data();
    end
    
    % Create main window
    fig = figure('Position', [100, 100, 1200, 700], ...
                'Name', 'EIDORS Interactive Simulation Interface', ...
                'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'CloseRequestFcn', @close_gui);
    
    % Create UI
    create_ui(fig);
    
    % Display startup information
    fprintf('EIDORS interactive interface startup completed\n');
    fprintf('=========================================\n');
    fprintf('Usage Instructions:\n');
    fprintf('1. Click on the left panel to place inclusions\n');
    fprintf('2. Adjust parameters: radius, conductivity, regularization parameter\n');
    fprintf('3. Click "Reconstruct" for standard reconstruction\n');
    fprintf('4. Click "High Quality Reconstruction" to test multiple EIDORS methods\n');
    fprintf('5. Use display buttons to view model mesh and true distribution\n');
    fprintf('=========================================\n');
    fprintf('Reconstruction Configuration (optimized by official code):\n');
    fprintf('- Forward model: %d elements (high resolution simulation)\n', size(app_data.fmdl_fwd.elems, 1));
    fprintf('- Reconstruction model: %d elements (coarse mesh reconstruction)\n', size(app_data.fmdl_rec.elems, 1));
    fprintf('- Default regularization parameter: %.0e (official recommended value)\n', 1e-1);
    fprintf('- Reconstruction method: inv_solve_diff_GN_one_step\n');
    fprintf('- Prior: prior_gaussian_HPF (official demo standard)\n');
    fprintf('=========================================\n');

    % Initialize application data (according to official code standards)
    function data = init_app_data()
        data = struct();
        data.inclusions = [];
        data.current_inclusion_idx = 0;
        
        try
            % n_elec = 16;
            % n_rings = 1;
            % options = {'no_meas_current','no_rotate_meas'};

            % % Forward model - high resolution (official recommended configuration)
            % params_fwd = mk_circ_tank(12, [], n_elec);
            % params_fwd.stimulation = mk_stim_patterns(n_elec, n_rings, '{ad}', '{ad}', options, 10);
            % params_fwd.solve = 'fwd_solve_1st_order';
            % params_fwd.system_mat = 'system_mat_1st_order';
            % data.fmdl_fwd = eidors_obj('fwd_model', params_fwd);
            % 
            % % Reconstruction model - coarse mesh
            % params_rec = mk_circ_tank(8, [], n_elec);
            % params_rec.stimulation = mk_stim_patterns(n_elec, n_rings, '{ad}', '{ad}', options, 10);
            % params_rec.solve = 'fwd_solve_1st_order';
            % params_rec.system_mat = 'system_mat_1st_order';
            % params_rec.jacobian = 'jacobian_adjoint';
            % data.fmdl_rec = eidors_obj('fwd_model', params_rec);
            
            % Use official template model
            imdl = mk_common_model('j2c2', 16);
            fmdl = imdl.fwd_model;
    
            data.fmdl_fwd = fmdl;
            data.fmdl_rec = imdl.fwd_model;
            data.inv_model_template = imdl;

            % background
            data.background_conductivity = 5;
    
            % Background image and simulation
            data.homg_img = mk_image(data.fmdl_fwd, data.background_conductivity);
            data.homg_data = fwd_solve(data.homg_img);
            
            % resize to normal
            data.uiR = 12;
    
            fprintf('EIDORS model initialization completed\n');
            fprintf('Forward model: %d elements, Reconstruction model: %d elements\n', ...
                    size(data.fmdl_fwd.elems, 1), size(data.fmdl_rec.elems, 1));
            fprintf('Uniform background conductivity: %.1f\n', data.homg_img.elem_data(1));
    
        catch ME
            fprintf('Initialization error: %s\n', ME.message);
        end
    end

    % Create user interface
    function create_ui(parent_fig)
        % Left panel
        panel_left = uipanel('Parent', parent_fig, ...
                            'Position', [0.02, 0.02, 0.48, 0.96], ...
                            'Title', 'Interactive Area - Click to place inclusions');
        
        app_data.ax_interact = axes('Parent', panel_left, ...
                                   'Position', [0.1, 0.1, 0.85, 0.85]);
        axis(app_data.ax_interact, 'equal');
        axis(app_data.ax_interact, [-15, 15, -15, 15]);
        title(app_data.ax_interact, 'Click to place circular inclusions');
        grid(app_data.ax_interact, 'on');
        
        set(app_data.ax_interact, 'ButtonDownFcn', @axes_click);
        set(app_data.ax_interact, 'HitTest', 'on');
        set(app_data.ax_interact, 'PickableParts', 'all');
        
        % Right upper panel - Control panel
        panel_control = uipanel('Parent', parent_fig, ...
                               'Position', [0.52, 0.52, 0.46, 0.46], ...
                               'Title', 'Control Panel');
        
        % Create controls
        create_controls(panel_control);
        
        % Right lower panel - Result display
        panel_result = uipanel('Parent', parent_fig, ...
                              'Position', [0.52, 0.02, 0.46, 0.48], ...
                              'Title', 'Reconstruction Results');
        
        app_data.ax_result = axes('Parent', panel_result, ...
                                 'Position', [0.1, 0.1, 0.85, 0.85]);
        axis(app_data.ax_result, 'equal');
        title(app_data.ax_result, 'Reconstruction results will be displayed here');
        
        % Initial display
        update_display();
    end

    % Create control panel
    function create_controls(parent_panel)
        % Radius control
        uicontrol('Parent', parent_panel, 'Style', 'text', ...
                 'Position', [10, 280, 100, 20], 'String', 'Inclusion Radius:');
        app_data.slider_radius = uicontrol('Parent', parent_panel, 'Style', 'slider', ...
                                          'Position', [120, 280, 200, 20], ...
                                          'Min', 0.5, 'Max', 4, 'Value', 2, ...
                                          'Callback', @update_radius_text);
        app_data.text_radius = uicontrol('Parent', parent_panel, 'Style', 'text', ...
                                        'Position', [330, 280, 50, 20], 'String', '2.0');
        
        % Conductivity control
        uicontrol('Parent', parent_panel, 'Style', 'text', ...
                 'Position', [10, 250, 100, 20], 'String', 'Conductivity Mult:');
        app_data.slider_conductivity = uicontrol('Parent', parent_panel, 'Style', 'slider', ...
                                                'Position', [120, 250, 200, 20], ...
                                                'Min', 0.1, 'Max', 10, 'Value', 5, ...
                                                'Callback', @update_conductivity_text);
        app_data.text_conductivity = uicontrol('Parent', parent_panel, 'Style', 'text', ...
                                              'Position', [330, 250, 50, 20], 'String', '5.0');
        
        % Regularization control (adjusted default value according to official code)
        uicontrol('Parent', parent_panel, 'Style', 'text', ...
                 'Position', [10, 220, 100, 20], 'String', 'Regularization:');
        app_data.slider_regularization = uicontrol('Parent', parent_panel, 'Style', 'slider', ...
                                                  'Position', [120, 220, 200, 20], ...
                                                  'Min', 1e-4, 'Max', 1e0, 'Value', 1e-1, ...
                                                  'Callback', @update_regularization_text);
        app_data.text_regularization = uicontrol('Parent', parent_panel, 'Style', 'text', ...
                                                'Position', [330, 220, 70, 20], 'String', '1.0e-1');
        
        % Main buttons
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [10, 180, 100, 30], 'String', 'Reconstruct', ...
                 'Callback', @perform_reconstruction_btn, ...
                 'BackgroundColor', [0.58, 0.65, 0.83]);
        
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [120, 180, 100, 30], 'String', 'Clear All', ...
                 'Callback', @clear_all_inclusions, ...
                 'BackgroundColor', [0.58, 0.65, 0.83]);
        
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [230, 180, 100, 30], 'String', 'Remove Last', ...
                 'Callback', @remove_last_inclusion, ...
                 'BackgroundColor', [0.58, 0.65, 0.83]);
        
        % Advanced buttons
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [10, 140, 120, 30], 'String', 'High Quality Recon', ...
                 'Callback', @advanced_reconstruction, ...
                 'BackgroundColor',[0.23, 0.35, 0.60]);
        
        app_data.checkbox_realtime = uicontrol('Parent', parent_panel, 'Style', 'checkbox', ...
                                              'Position', [140, 140, 180, 20], ...
                                              'String', 'Real-time Recon (slower)', 'Value', 0);
        
        % Display button group
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [10, 100, 100, 25], 'String', 'Show Forward Model', ...
                 'Callback', @show_forward_model, ...
                 'BackgroundColor', [0.5, 0.5, 0.8]);
        
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [120, 100, 100, 25], 'String', 'Show Recon Model', ...
                 'Callback', @show_reconstruction_model, ...
                 'BackgroundColor', [0.5, 0.5, 0.8]);
        
        uicontrol('Parent', parent_panel, 'Style', 'pushbutton', ...
                 'Position', [230, 100, 100, 25], 'String', 'Show True Distrib', ...
                 'Callback', @show_true_distribution, ...
                 'BackgroundColor', [0.5, 0.5, 0.8]);
        
        
        % Information display
        app_data.text_info = uicontrol(...
            'Parent', parent_panel, ...
            'Style', 'listbox', ...
            'Position', [10, 10, 370, 80], ...
            'Max', 2, 'Min', 0, ...
            'HorizontalAlignment', 'left', ...
            'Enable', 'inactive', ...
            'BackgroundColor', 'white' ...
        );
    end

    % Click callback
    function axes_click(~, ~)
        try
            point = get(app_data.ax_interact, 'CurrentPoint');
            x = point(1, 1);
            y = point(1, 2);
            
            fprintf('Click position: x=%.2f, y=%.2f\n', x, y);
            
            distance_from_center = sqrt(x^2 + y^2);
            if distance_from_center > 11
                msg = sprintf('Click position outside boundary (distance from center: %.2f)', distance_from_center);
                fprintf('%s\n', msg);
                update_info(msg);
                return;
            end
            
            app_data.current_inclusion_idx = app_data.current_inclusion_idx + 1;
            inclusion = struct();
            inclusion.x = x;
            inclusion.y = y;
            inclusion.radius = get(app_data.slider_radius, 'Value');
            inclusion.conductivity = get(app_data.slider_conductivity, 'Value');
            inclusion.id = app_data.current_inclusion_idx;
            
            app_data.inclusions = [app_data.inclusions, inclusion];
            
            msg = sprintf('Added inclusion #%d at position (%.2f, %.2f)', inclusion.id, x, y);
            fprintf('%s\n', msg);
            update_info(msg);
            
            update_display();
            
            if get(app_data.checkbox_realtime, 'Value')
                perform_reconstruction();
            end
            
        catch ME
            fprintf('Click handling error: %s\n', ME.message);
        end
    end

    % Update display
    function update_display()
        try
            cla(app_data.ax_interact);
            
            set(app_data.ax_interact, 'ButtonDownFcn', @axes_click);
            set(app_data.ax_interact, 'HitTest', 'on');
            set(app_data.ax_interact, 'PickableParts', 'all');
            
            draw_tank_boundary();
            
            for i = 1:length(app_data.inclusions)
                inc = app_data.inclusions(i);
                draw_inclusion(inc);
            end
            
            axis(app_data.ax_interact, [-15, 15, -15, 15]);
            axis(app_data.ax_interact, 'equal');
            grid(app_data.ax_interact, 'on');
            title(app_data.ax_interact, sprintf('Inclusion count: %d (click to add)', length(app_data.inclusions)));
            
        catch ME
            fprintf('Display update error: %s\n', ME.message);
        end
    end

    % Draw boundary
    function draw_tank_boundary()
        theta = linspace(0, 2*pi, 100);
        x_boundary = 12 * cos(theta);
        y_boundary = 12 * sin(theta);
        h_boundary = plot(app_data.ax_interact, x_boundary, y_boundary, 'k-', 'LineWidth', 2);
        hold(app_data.ax_interact, 'on');
        set(h_boundary, 'HitTest', 'off');
        
        n_elec = 16;
        for i = 1:n_elec
            angle = 2*pi*(i-1)/n_elec;
            x_elec = 12 * cos(angle);
            y_elec = 12 * sin(angle);
            h_elec = plot(app_data.ax_interact, x_elec, y_elec, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            h_text = text(app_data.ax_interact, x_elec*1.1, y_elec*1.1, num2str(i), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
            set(h_elec, 'HitTest', 'off');
            set(h_text, 'HitTest', 'off');
        end
    end

    % Draw inclusion
    function draw_inclusion(inc)
        theta = linspace(0, 2*pi, 50);
        x_inc = inc.x + inc.radius * cos(theta);
        y_inc = inc.y + inc.radius * sin(theta);

        bg_cond = app_data.homg_img.elem_data(1);
        
        if inc.conductivity > bg_cond
            color = 'r';
        else
            color = 'b';
        end
        
        h_fill = fill(app_data.ax_interact, x_inc, y_inc, color, 'FaceAlpha', 0.5, 'EdgeColor', color);
        h_text = text(app_data.ax_interact, inc.x, inc.y, num2str(inc.id), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white');
        set(h_fill, 'HitTest', 'off');
        set(h_text, 'HitTest', 'off');
    end

    % Perform reconstruction (enhanced debug version)
    function perform_reconstruction()
        if isempty(app_data.inclusions)
            msg = 'No inclusions, cannot reconstruct';
            fprintf('%s\n', msg);
            update_info(msg);
            return;
        end
        
        try
            fprintf('=== Starting reconstruction process ===\n');
            update_info('Performing reconstruction...');
            
            % Create inclusion image
            inh_img = create_inclusion_image();
            
            % Forward solve
            fprintf('Forward solving inclusion data...\n');
            inh_data = fwd_solve(inh_img);
            
            % Calculate data difference
            data_diff = inh_data.meas - app_data.homg_data.meas;
            signal_strength = norm(data_diff);
            max_diff = max(abs(data_diff));
            mean_diff = mean(abs(data_diff));
            
            fprintf('Signal analysis:\n');
            fprintf('  Signal strength (L2 norm): %.6f\n', signal_strength);
            fprintf('  Maximum difference: %.6f\n', max_diff);
            fprintf('  Mean difference: %.6f\n', mean_diff);
            fprintf('  Number of measurements: %d\n', length(data_diff));
            
            % Check if signal is strong enough
            if signal_strength < 1e-10
                warning('Signal strength too small, may not be able to reconstruct');
            end
            
            % Create inverse problem model
            inv_mdl = create_inverse_model();
            
            % Perform reconstruction (note parameter order: model, inclusion data, uniform data)
            fprintf('Performing inverse reconstruction...\n');
            img_rec = inv_solve(inv_mdl, inh_data, app_data.homg_data);
            
            if isempty(img_rec) || ~isfield(img_rec, 'elem_data')
                error('Reconstruction failed: unable to generate valid reconstruction image');
            end
            
            img_rec.elem_data = - img_rec.elem_data;
            img_rec.elem_data = img_rec.elem_data + app_data.homg_img.elem_data;
            
            % Analyze reconstruction results
            rec_data = img_rec.elem_data;
            fprintf('Reconstruction result analysis:\n');
            fprintf('  Number of reconstruction elements: %d\n', length(rec_data));
            fprintf('  Value range: [%.6f, %.6f]\n', min(rec_data), max(rec_data));
            fprintf('  Standard deviation: %.6f\n', std(rec_data));
            fprintf('  Non-zero elements: %d\n', sum(abs(rec_data) > 1e-10));
            
            % Display results
            display_reconstruction_result(img_rec);
            show_reconstruction_stats(img_rec);
            
            fprintf('=== Reconstruction completed ===\n');
            update_info(sprintf('Reconstruction completed - Signal strength:%.2e, Recon range:[%.3f,%.3f]', ...
                               signal_strength, min(rec_data), max(rec_data)));
            
        catch ME
            fprintf('Reconstruction error: %s\n', ME.message);
            fprintf('Error stack:\n');
            for i = 1:length(ME.stack)
                fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
            end
            update_info(['Reconstruction error: ' ME.message]);
        end
    end

    % Create inclusion image (detailed debug version)
    function inh_img = create_inclusion_image()
        % First create image with uniform background conductivity
        inh_img = mk_image(app_data.fmdl_fwd, app_data.background_conductivity);
        if isempty(app_data.inclusions)
            return;
        end

        % Calculate element centers
        elem_centers = calc_elem_centers();
        elem_count   = size(app_data.fmdl_fwd.elems, 1);

        % —— Calculate the ratio between UI coordinate system and model coordinate system —— 
        uiR    = app_data.uiR;  % Container radius drawn in UI
        % Actual maximum radius of model (calculated from all nodes)
        node_xy = app_data.fmdl_fwd.nodes(:,1:2);
        modelR  = max( sqrt(sum(node_xy.^2, 2)) );
        % ————————————————————————————————————————————————

        fprintf('=== Inclusion setup details ===\n');
        total_affected = 0;

        for i = 1:length(app_data.inclusions)
            inc = app_data.inclusions(i);
            fprintf('Processing inclusion #%d: UI position(%.2f,%.2f), UI radius%.2f, conductivity%.2f\n', ...
                    inc.id, inc.x, inc.y, inc.radius, inc.conductivity);

            % —— Scale UI units (x,y,radius) to model units —— 
            x0 = inc.x      / uiR * modelR;
            y0 = inc.y      / uiR * modelR;
            r0 = inc.radius / uiR * modelR;
            fprintf('  Mapped model coordinates: (%.3f,%.3f), radius%.3f\n', x0, y0, r0);
            % ——————————————————————————————————————————————

            % Calculate distance from each element center to inclusion center
            distances = sqrt((elem_centers(:,1)-x0).^2 + (elem_centers(:,2)-y0).^2);

            % Strict/loose mapping
            inside_strict = distances <= r0;
            inside_loose  = distances <= r0 * 1.5;

            % If strict condition doesn't match, use loose condition
            if ~any(inside_strict)
                fprintf('  Elements within strict radius: 0, using loose condition\n');
                inside_elements = inside_loose;
            else
                fprintf('  Elements within strict radius: %d\n', sum(inside_strict));
                inside_elements = inside_strict;
            end

            % If still no match, force select nearest 5 elements
            if ~any(inside_elements)
                fprintf('  Loose condition also no elements, forcing selection of nearest 5\n');
                [~, idx_sorted] = sort(distances);
                force_n = min(5, elem_count);
                inside_elements = false(elem_count, 1);
                inside_elements(idx_sorted(1:force_n)) = true;
            end

            % Find finally affected element indices
            affected_indices = find(inside_elements);
            n_aff = length(affected_indices);
            total_affected = total_affected + n_aff;

            % Apply conductivity
            inh_img.elem_data(affected_indices) = inc.conductivity;

            fprintf('  Finally affected elements: %d (first 10 indices: %s)\n', ...
                    n_aff, mat2str(affected_indices(1:min(10,n_aff))));
            fprintf('  Set conductivity: %.2f\n', inc.conductivity);
        end

        fprintf('Total affected elements: %d / %d\n', total_affected, elem_count);

        % Verify and print conductivity distribution
        unique_vals = unique(inh_img.elem_data);
        fprintf('Conductivity distribution values: %s\n', mat2str(unique_vals, 3));
        for v = unique_vals'
            cnt = sum(inh_img.elem_data == v);
            fprintf('  Conductivity %.2f: %d elements\n', v, cnt);
        end
        fprintf('==================\n');
    end

    % Calculate element centers (fixed version)
    function centers = calc_elem_centers()
        nodes = app_data.fmdl_fwd.nodes;
        elems = app_data.fmdl_fwd.elems;
        n_elems = size(elems, 1);
        
        % Check node dimensions
        if size(nodes, 2) >= 2
            node_coords = nodes(:, 1:2);  % Take only x,y coordinates
        else
            error('Insufficient node coordinate dimensions');
        end
        
        centers = zeros(n_elems, 2);
        
        % Calculate center point of each element
        for i = 1:n_elems
            elem_nodes = elems(i, :);
            % Remove possible 0 nodes (some mesh formats)
            elem_nodes = elem_nodes(elem_nodes > 0);
            if length(elem_nodes) >= 3
                centers(i, :) = mean(node_coords(elem_nodes, :), 1);
            end
        end
        
        % Debug information: display mesh range
        fprintf('Mesh information:\n');
        fprintf('  Number of nodes: %d, Number of elements: %d\n', size(nodes,1), n_elems);
        fprintf('  Node coordinate range: x[%.2f, %.2f], y[%.2f, %.2f]\n', ...
               min(node_coords(:,1)), max(node_coords(:,1)), ...
               min(node_coords(:,2)), max(node_coords(:,2)));
        fprintf('  Element center range: x[%.2f, %.2f], y[%.2f, %.2f]\n', ...
               min(centers(:,1)), max(centers(:,1)), ...
               min(centers(:,2)), max(centers(:,2)));
        
        % Display sample element centers
        sample_indices = [1, round(n_elems/4), round(n_elems/2), round(3*n_elems/4), n_elems];
        fprintf('  Sample element centers: ');
        for idx = sample_indices
            if idx <= n_elems
                fprintf('[%.2f,%.2f] ', centers(idx,1), centers(idx,2));
            end
        end
        fprintf('\n');
    end

    % Create inverse problem model (optimized according to official code)
    function inv_mdl = create_inverse_model()
        % Optimized configuration according to official tutorial
        clear inv2d;
        inv2d.name = 'EIT inverse';
        inv2d.solve = 'inv_solve_diff_GN_one_step';  % Official standard solver
        inv2d.hyperparameter.value = get(app_data.slider_regularization, 'Value');
        inv2d.RtR_prior = 'prior_gaussian_HPF';      % Official recommended high-pass filter prior
        inv2d.reconst_type = 'difference';
        inv2d.jacobian_bkgnd.value = 1;
        inv2d.fwd_model = app_data.fmdl_rec;
        
        inv_mdl = eidors_obj('inv_model', inv2d);
        
        fprintf('Reconstruction configuration: solver=%s, prior=%s, regularization=%.0e\n', ...
               inv2d.solve, inv2d.RtR_prior, inv2d.hyperparameter.value);
    end

    % Display reconstruction results
    function display_reconstruction_result(img_rec)
        if ~isfield(img_rec, 'elem_data') || isempty(img_rec.elem_data)
            fprintf('Reconstruction result has no data\n');
            return;
        end
        
        axes(app_data.ax_result);
        cla(app_data.ax_result);
        
        try
            show_slices(img_rec);
            title(app_data.ax_result, 'Reconstructed conductivity change');
            
            rec_data = img_rec.elem_data;
            data_range = [min(rec_data), max(rec_data)];
            
            if abs(data_range(2) - data_range(1)) < 1e-10
                max_abs = max(abs(rec_data));
                if max_abs > 0
                    data_range = [-max_abs, max_abs];
                end
            end
            
            if data_range(1) ~= data_range(2)
                set(app_data.ax_result, 'CLim', data_range);
            end
            
            colorbar(app_data.ax_result);
            colormap(app_data.ax_result, 'jet');
            
            fprintf('Reconstruction result display: range[%.6f, %.6f], std=%.6f\n', ...
                   data_range(1), data_range(2), std(rec_data));
            
        catch ME1
            fprintf('show_slices failed: %s\n', ME1.message);
            try
                display_with_patch(img_rec);
            catch ME2
                fprintf('All display methods failed: %s\n', ME2.message);
            end
        end
    end

    % Display using patch
    function display_with_patch(img_rec)
        nodes = img_rec.fwd_model.nodes;
        elems = img_rec.fwd_model.elems;
        elem_data = img_rec.elem_data;
        
        if size(nodes, 2) > 2
            nodes = nodes(:, 1:2);
        end
        
        patch('Vertices', nodes, 'Faces', elems, ...
              'FaceVertexCData', elem_data, 'FaceColor', 'flat', ...
              'EdgeColor', 'none', 'Parent', app_data.ax_result);
        
        axis(app_data.ax_result, 'equal');
        title(app_data.ax_result, 'Reconstructed conductivity change');
        colorbar(app_data.ax_result);
        colormap(app_data.ax_result, 'jet');
        
        max_val = max(abs(elem_data));
        if max_val > 0
            set(app_data.ax_result, 'CLim', [-max_val, max_val]);
        end
    end

    % Display statistics
    function show_reconstruction_stats(img_rec)
        try
            rec_data = img_rec.elem_data;
            max_rec = max(abs(rec_data));
            min_rec = min(rec_data);
            mean_rec = mean(rec_data);
            std_rec = std(rec_data);
            
            fprintf('=== Reconstruction Statistics ===\n');
            fprintf('Number of data points: %d\n', length(rec_data));
            fprintf('Maximum value: %.4f\n', max_rec);
            fprintf('Minimum value: %.4f\n', min_rec);
            fprintf('Mean value: %.4f\n', mean_rec);
            fprintf('Standard deviation: %.4f\n', std_rec);
            fprintf('Dynamic range: %.4f\n', max_rec - min_rec);
            fprintf('================\n');
            
        catch ME
            fprintf('Statistics calculation error: %s\n', ME.message);
        end
    end

    % Update information display
    function update_info(msg)
        try
 
            t = char(datetime('now','Format','HH:mm:ss'));
            new_entry = sprintf('%s - %s', t, msg);

            old = get(app_data.text_info, 'String');
            if ischar(old)
                old = cellstr(old);
            end
    
            old{end+1} = new_entry;
            max_lines = 100;
            if numel(old) > max_lines
                old = old(end-max_lines+1:end);
            end
    
            set(app_data.text_info, ...
                'String', old, ...
                'Value', numel(old), ...     
                'ListboxTop', numel(old));    
    
            drawnow; 
    
        catch
            fprintf('%s - %s\n', char(datetime('now','Format','HH:mm:ss')), msg);
        end
    end

    % UI callback functions
    function update_radius_text(~, ~)
        val = get(app_data.slider_radius, 'Value');
        set(app_data.text_radius, 'String', sprintf('%.1f', val));
    end

    function update_conductivity_text(~, ~)
        val = get(app_data.slider_conductivity, 'Value');
        set(app_data.text_conductivity, 'String', sprintf('%.1f', val));
    end

    function update_regularization_text(~, ~)
        val = get(app_data.slider_regularization, 'Value');
        % Use scientific notation display, consistent with official code
        set(app_data.text_regularization, 'String', sprintf('%.1e', val));
    end

    function perform_reconstruction_btn(~, ~)
        perform_reconstruction();
    end

    function clear_all_inclusions(~, ~)
        app_data.inclusions = [];
        app_data.current_inclusion_idx = 0;
        update_display();
        cla(app_data.ax_result);
        title(app_data.ax_result, 'Reconstruction results will be displayed here');
        fprintf('All inclusions cleared\n');
        update_info('All inclusions cleared');
    end

    function remove_last_inclusion(~, ~)
        if ~isempty(app_data.inclusions)
            app_data.inclusions(end) = [];
            update_display();
            fprintf('Last inclusion removed, %d remaining\n', length(app_data.inclusions));
            update_info(sprintf('Last inclusion removed, %d remaining', length(app_data.inclusions)));
            
            if get(app_data.checkbox_realtime, 'Value') && ~isempty(app_data.inclusions)
                perform_reconstruction();
            elseif isempty(app_data.inclusions)
                cla(app_data.ax_result);
                title(app_data.ax_result, 'Reconstruction results will be displayed here');
            end
        else
            fprintf('No inclusions to remove\n');
            update_info('No inclusions to remove');
        end
    end

    function advanced_reconstruction(~, ~)
        if isempty(app_data.inclusions)
            fprintf('No inclusions, cannot reconstruct\n');
            update_info('No inclusions, cannot reconstruct');
            return;
        end
        
        try
            fprintf('=== High Quality Reconstruction: Testing Multiple Official Standard Methods ===\n');
            update_info('Performing high quality reconstruction...');
            
            inh_img = create_inclusion_image();
            inh_data = fwd_solve(inh_img);
            
            % Configure multiple methods according to official code
            methods = struct();
            
            % Method 1: Gaussian HPF prior (used in official demo)
            methods(1).name = 'Gaussian_HPF';
            methods(1).solve = 'inv_solve_diff_GN_one_step';
            methods(1).RtR_prior = 'prior_gaussian_HPF';
            methods(1).hyperparameter = 1e-1;
            
            % Method 2: TV prior
            methods(2).name = 'TV_prior';
            methods(2).solve = 'inv_solve_diff_GN_one_step';
            methods(2).RtR_prior = 'prior_TV';
            methods(2).hyperparameter = 1e-2;
            
            % Method 3: Laplace prior
            methods(3).name = 'Laplace';
            methods(3).solve = 'inv_solve_diff_GN_one_step';
            methods(3).RtR_prior = 'prior_laplace';
            methods(3).hyperparameter = 1e-3;
            
            % Method 4: NOSER prior
            methods(4).name = 'NOSER';
            methods(4).solve = 'inv_solve_diff_GN_one_step';
            methods(4).RtR_prior = 'prior_noser';
            methods(4).hyperparameter = 1e-2;
            
            best_img = [];
            best_score = -inf;
            best_method = '';
            
            for i = 1:length(methods)
                try
                    fprintf('Testing method %d/%d: %s\n', i, length(methods), methods(i).name);
                    
                    % Create inverse problem model
                    clear inv2d;
                    inv2d.name = sprintf('EIT_%s', methods(i).name);
                    inv2d.solve = methods(i).solve;
                    inv2d.hyperparameter.value = methods(i).hyperparameter;
                    inv2d.RtR_prior = methods(i).RtR_prior;
                    inv2d.reconst_type = 'difference';
                    inv2d.jacobian_bkgnd.value = 1;
                    inv2d.fwd_model = app_data.fmdl_rec;
                    
                    inv_mdl = eidors_obj('inv_model', inv2d);
                    
                    % Perform reconstruction
                    img_rec = inv_solve(inv_mdl, inh_data, app_data.homg_data);
                    
                    if ~isempty(img_rec) && isfield(img_rec, 'elem_data')
                        img_rec.elem_data = - img_rec.elem_data;
                        img_rec.elem_data = img_rec.elem_data + app_data.homg_img.elem_data;
                        
                        rec_data = img_rec.elem_data;
                        
                        % Use more reasonable scoring criteria
                        % 1. Check dynamic range
                        dynamic_range = max(rec_data) - min(rec_data);
                        % 2. Check standard deviation
                        std_score = std(rec_data);
                        % 3. Check non-zero element ratio
                        nonzero_ratio = sum(abs(rec_data) > 1e-6) / length(rec_data);
                        
                        % Composite score
                        score = dynamic_range * std_score * (1 + nonzero_ratio);
                        
                        fprintf('  Dynamic range: %.6f, Std dev: %.6f, Non-zero ratio: %.3f\n', ...
                               dynamic_range, std_score, nonzero_ratio);
                        fprintf('  Composite score: %.6f\n', score);
                        
                        if score > best_score && ~any(isnan(rec_data)) && ~any(isinf(rec_data))
                            best_score = score;
                            best_img = img_rec;
                            best_method = methods(i).name;
                        end
                    end
                    
                catch ME
                    fprintf('  Method %s failed: %s\n', methods(i).name, ME.message);
                end
            end
            
            if ~isempty(best_img)
                display_reconstruction_result(best_img);
                show_reconstruction_stats(best_img);
                fprintf('=== High Quality Reconstruction Completed ===\n');
                fprintf('Best method: %s, Score: %.6f\n', best_method, best_score);
                update_info(sprintf('Best method: %s, Score: %.6f', best_method, best_score));
            else
                fprintf('All reconstruction methods failed\n');
                update_info('All methods failed');
            end
            
        catch ME
            fprintf('High quality reconstruction error: %s\n', ME.message);
            update_info(['Reconstruction error: ' ME.message]);
        end
    end

    function show_forward_model(~, ~)
        try
            figure('Name', 'Forward Model Mesh', 'NumberTitle', 'off');
            show_fem(app_data.fmdl_fwd);
            title('Forward Model Mesh (High Resolution)');
            fprintf('Forward model: %d elements, %d nodes\n', ...
                   size(app_data.fmdl_fwd.elems, 1), size(app_data.fmdl_fwd.nodes, 1));
        catch ME
            fprintf('Forward model display error: %s\n', ME.message);
        end
    end

    function show_reconstruction_model(~, ~)
        try
            figure('Name', 'Reconstruction Model Mesh', 'NumberTitle', 'off');
            show_fem(app_data.fmdl_rec);
            title('Reconstruction Model Mesh (Coarse Mesh)');
            fprintf('Reconstruction model: %d elements, %d nodes\n', ...
                   size(app_data.fmdl_rec.elems, 1), size(app_data.fmdl_rec.nodes, 1));
        catch ME
            fprintf('Reconstruction model display error: %s\n', ME.message);
        end
    end

    function show_true_distribution(~, ~)
        if isempty(app_data.inclusions)
            fprintf('No inclusions to display\n');
            update_info('No inclusions to display');
            return;
        end
        
        try
            % Create inclusion image
            true_img = create_inclusion_image();
            
            % Display true distribution in result area
            axes(app_data.ax_result);
            cla(app_data.ax_result);
            
            % Check if data has variation
            unique_vals = unique(true_img.elem_data);
            fprintf('True distribution conductivity values: %s\n', mat2str(unique_vals, 3));
            
            if isscalar(unique_vals)
                % If only one value, inclusion setting failed
                text(app_data.ax_result, 0.5, 0.5, {'Inclusion mapping failed!', 'Check position and size', 'or click "Debug Mesh"'}, ...
                    'HorizontalAlignment', 'center', 'Units', 'normalized', ...
                    'FontSize', 12, 'Color', 'red');
                fprintf('Error: Inclusions failed to map to any mesh elements\n');
                update_info('Error: Inclusion mapping failed');
                return;
            end
            
            % Display true distribution
            show_slices(true_img);
            title(app_data.ax_result, 'True Conductivity Distribution');
            colorbar(app_data.ax_result);
            colormap(app_data.ax_result, 'jet');
            
            % Set appropriate color range
            clim_val = [min(true_img.elem_data), max(true_img.elem_data)];
            if clim_val(1) ~= clim_val(2)
                set(app_data.ax_result, 'CLim', clim_val);
            end
            
            fprintf('True conductivity distribution display completed\n');
            update_info('True distribution displayed');
            
        catch ME
            fprintf('True distribution display error: %s\n', ME.message);
            update_info(['Display error: ' ME.message]);
        end
    end

    function close_gui(~, ~)
        selection = questdlg('Are you sure you want to close the EIDORS interactive interface?', ...
                           'Confirm Close', 'Yes', 'No', 'No');
        if strcmp(selection, 'Yes')
            clear app_data;
            delete(gcf);
        end
    end

end