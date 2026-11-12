function [CenXY_plot, Fall, R] = ProlongedFieldAbsenceTest(path, path2, start_frame, full, threshold, timeR, step)
    % ProlongedFieldAbsenceTest Analyze the prolonged field absence test with eight houses.
    %   [CenXY_plot, Fall, R] = ProlongedFieldAbsenceTest(path, path2, start_frame, ...
    %       full, threshold, timeR, step)
    %
    % Inputs
    %   path        Absolute folder containing 1024 x 1280 grayscale .tiff images.
    %               Filenames must end in an underscore followed by a frame number.
    %   path2       Absolute output folder for cached data, results, and figures.
    %   start_frame First image to analyze in the sorted file list (1 based).
    %   full        Recording duration in hours; 0 uses all remaining images [0].
    %   threshold   Maximum displacement per analyzed interval, in pixels [10].
    %   timeR       Time between source images, in seconds [10].
    %   step        Analyze every step-th source image [1]. The analyzed interval
    %               is step * timeR seconds. Set this to match each recording.
    %
    % Fixed setup
    %   Image size is [1024, 1280] = [height, width]; all 8 houses are analyzed.
    %   Stationary bouts require 450 seconds at or below the movement threshold.
    %   Spatial calibration is 100 pixels per 7 cm. Defaults give 7 mm/interval.
    %
    % Outputs
    %   CenXY_plot  Analyzed positions in pixels, with the vertical axis flipped.
    %   Fall        Binary tracking masks retained every 10 analyzed frames.
    %   R           Structure containing displacement, stationary occupancy,
    %               eight-house preference, and directional/axial vector results.
    %               R.PI: all-frame occupancy fractions among the eight houses.
    %               R.In_house_ratio: fraction of all frames inside any house.
    %               R.RIH: stationary house frames / all analyzed frames.
    %               R.RRE: stationary frames / all analyzed frames.
    %               R.pd: direction in radians; R.po: doubled-angle axial result.
    %
    % Requirements: MATLAB Image Processing Toolbox, rgb.m, and eight house ROIs
    % with a Position polygon in each Houses{j,1}. Edit house_roi_file below if needed.
    % Cached MAT files are reused when present. Use a new output folder when
    % changing the source images, frame selection, or house polygons.
    % The function changes the working directory to path2 and closes figures.
    % Author: Li Zhang. Original version: 2021-11-22.

    %% 1. Fixed setup and optional arguments
    resolution = [1024, 1280]; % [height, width] in pixels
    house_roi_file = 'Z:\MatlabFunctions\Behavior\Avg_Houses.mat';
    if nargin < 4
        full = 0;
    end
    if nargin < 5
        threshold = 10;
    end
    if nargin < 6
        timeR = 10;
    end
    if nargin < 7
        step = 1;
    end

    % Number of consecutive displacement samples required for a stationary bout.
    fstd = 450 / (step * timeR); % 450sec (7.5min)

    %% 2. Output folder and recording label
    if ~exist(path2, 'dir')
        mkdir(path2);
    end
    cd(path2);

    PathSlash = find(path == '\', 1, 'last');
    Mid0 = path(PathSlash + 1:end);
    Mid = Mid0;
    Mid(Mid0 == '_') = ' ';

    %% 3. Load or collect image paths in frame-number order
    if exist(fullfile(path2, 'All_frames_path.mat'), 'file')  == 2
        load(fullfile(path2, 'All_frames_path.mat'), 'x');
        disp('Loading cached image paths...');
    else
        x = dir(fullfile(path, '*.tiff'));
        order = zeros(length(x), 1);
        for i = 1:length(x)
            a1 = find(x(i).name == '_', 1, 'last');
            a2 = find(x(i).name == '.', 1, 'first');
            order(i, 1) = str2double(x(i).name(a1 + 1:a2 - 1));
        end
        [~, order] = sort(order);
        order = order - min(order) + 1;
        x = x(order);
        save(fullfile(path2, 'All_frames_path.mat'), 'x');
    end

    x = x(start_frame:step:end);
    fprintf('Analyzing every %d source images (%g seconds per interval).\n', step, step * timeR);

    %% 4. Estimate the background from every tenth image, up to analyzed image 1800
    if exist(fullfile(path2, 'Averaged_F.mat'), 'file')  == 2
        load(fullfile(path2, 'Averaged_F.mat'), 'F_averaged');
        disp('Loading cached background image...');
    else
        if length(x) < 1800
            f_id = 1:10:length(x);
        else
            f_id = 1:10:1800;
        end
        f_id_all = zeros(length(f_id), resolution(1), resolution(2));
        for f_id_num = 1:length(f_id)
            t = Tiff([x(f_id(f_id_num)).folder, '\', x(f_id(f_id_num)).name]);
            f_id_all(f_id_num, :, :) = read(t);
        end
        F_averaged = squeeze(mean(f_id_all, 1));
        save(fullfile(path2, 'Averaged_F.mat'), 'F_averaged');
    end

    %% 5. Load the eight house polygons
    if exist(fullfile(path2, 'Houses_ROI.mat'), 'file') == 2
        load(fullfile(path2, 'Houses_ROI.mat'), 'Houses');
        disp('Loading cached house polygons...');
    else
        figure('position', [2 161 1052 835]);
        imagesc(F_averaged);
        colormap(gray);
        load(house_roi_file, 'Avg_Houses');
        Houses = Avg_Houses;
        save(fullfile(path2, 'Houses_ROI.mat'), 'Houses');
    end

    %% 6. Track the mouse and retain sparse masks for occupancy heatmaps
    if exist(fullfile(path2, 'CenXY_plot.mat'), 'file')  == 2
        load(fullfile(path2, 'CenXY_plot.mat'), 'CenXY_plot', 'Fall');
        disp('Loading cached tracking results...');
    else
        t_start_tracking = tic;
        Fall = zeros(floor((length(x) - 1) / 10), resolution(1), resolution(2)); % one mask per 10 analyzed frames
        Fall = uint8(Fall);
        Fall_c = 0;
        CenXY = zeros(length(x), 2);
        for i = 1:length(x)
            t = Tiff([x(i).folder, '\', x(i).name]);
            F1 = read(t);
            % Subtract the frame from the background, then remove small structures.
            F_averaged = uint8(F_averaged);
            F0 = F_averaged - F1;
            F0(F0 < mean(mean(F0)) * 10) = 0;
            F0(F0 > 0) = 1;
            se = strel('disk', 10);
            F01 = imopen(F0, se);
            se = strel('disk', 20);
            F02 = imclose(F01, se);

            [L, numComponents] = bwlabel(F02);
            % Keep the largest connected component when several are detected.
            if numComponents > 1
                disp(['multiple components detected, # ', num2str(i), ' has ', num2str(numComponents), ' components']);

                componentStats = regionprops(L, 'Area');
                allAreas = [componentStats.Area];
                [~, largestComponentIdx] = max(allAreas);
                largestComponent = (L == largestComponentIdx);

                F02 = uint8(largestComponent);
                disp('The largest component has been kept...');

            end

            % First recovery: fill holes and repeat the standard morphology.
            if numComponents == 0
                F0 = imfill(F0, 'holes');
                se = strel('disk', 10);
                F01 = imopen(F0, se);
                se = strel('disk', 20);
                F02 = imclose(F01, se);
                [~, numComponents] = bwlabel(F02);
            end

            % Second recovery: use a smaller opening and save diagnostic images.
            if numComponents == 0
                disp('No components left!!!');
                if ~exist(fullfile(path2, 'no components'), 'dir')
                    mkdir(fullfile(path2, 'no components'));
                end
                h_check = figure('position', [327 153 750 750], 'visible', 'off');
                subplot(2, 2, 1);
                imagesc(F1);
                subplot(2, 2, 2);
                imagesc(F0);
                subplot(2, 2, 3);
                imagesc(F01);
                subplot(2, 2, 4);
                imagesc(F02);
                hold on;
                plot(CenXY(i, 1), CenXY(i, 2), 'or');
                colormap('gray');
                saveas(h_check, fullfile(path2, ['no components\', 'frame_', num2str(i), '.png']));
                close(h_check);
                se = strel('disk', 1); % one-pixel opening in the recovery pass
                F01 = imopen(F0, se);
                se = strel('disk', 20);
                F02 = imclose(F01, se);
                [L, numComponents] = bwlabel(F02);
                if numComponents > 1
                    disp(['multiple components detected, # ', num2str(i), ' has ', num2str(numComponents), ' components']);

                    componentStats = regionprops(L, 'Area');
                    allAreas = [componentStats.Area];
                    [~, largestComponentIdx] = max(allAreas);
                    largestComponent = (L == largestComponentIdx);

                    F02 = uint8(largestComponent);
                    disp('The largest component has been kept...');
                end
                h_check = figure('position', [327 153 750 750], 'visible', 'on');
                subplot(2, 2, 1);
                imagesc(F1);
                subplot(2, 2, 2);
                imagesc(F0);
                subplot(2, 2, 3);
                imagesc(F01);
                subplot(2, 2, 4);
                imagesc(F02);
                hold on;
                plot(CenXY(i, 1), CenXY(i, 2), 'or');
                colormap('gray');
                saveas(h_check, fullfile(path2, ['no components\', 'frame_', num2str(i), '_fixed.png']));
                close(h_check);
            end

            % If detection still fails, retain the previous centroid.
            stat = regionprops(F02);
            if isempty(stat)
                disp(['Frame: ', num2str(i), ' has no detected component.']);
                CenXY(i, 1) = CenXY(i - 1, 1);
                CenXY(i, 2) = CenXY(i - 1, 2);
            else
                CenXY(i, 1) = stat.Centroid(1);
                CenXY(i, 2) = stat.Centroid(2);
            end
            if mod(i, 10) == 0
                Fall_c = Fall_c + 1;
                Fall(Fall_c, :, :) = F02;
                if mod(i, 300) == 0
                    fprintf('Tracking mouse positions: frame %d\n', i);
                end
            end
        end
        CenXY_plot(:, 1) = CenXY(:, 1);
        CenXY_plot(:, 2) = resolution(1) - CenXY(:, 2);
        save(fullfile(path2, 'CenXY_plot.mat'), 'CenXY', 'CenXY_plot', 'Fall');
        t_end_tracking = toc(t_start_tracking);
        fprintf('Tracking mouse positions took: %f seconds\n', t_end_tracking);
    end
    %% 7. Apply the optional duration limit
    if full ~= 0 % optionally restrict the analyzed duration
        max_length = full * 3600 / (step * timeR); % how many datapoints during n hrs
        if length(CenXY_plot) >= max_length
            CenXY_plot = CenXY_plot(1:max_length, :);
            Fall = Fall(1:max_length / 10, :, :);
        else
            disp(['full length is less than ', full, ' hrs']);
        end
    end
    %% 8. Identify sustained stationary periods
    % Retain the original sample indexing and 450-second moving-window rule.
    diffs = diff(CenXY_plot, 1, 1);        % Differences between consecutive rows
    dis_trace = sqrt(sum(diffs.^2, 2));   % Euclidean distance for each row
    id_state = dis_trace <= threshold;

    id_fstd = movmean(id_state, fstd, 'Endpoints', 'discard');
    id_statef = find(id_fstd == 1);
    state_id = [];
    for i = 1:fstd
        state_id = [state_id; id_statef + i - 1];
    end
    state_id = unique(state_id);
    all_id = (1:length(CenXY_plot))';

    state_id(state_id > size(CenXY_plot, 1)) = [];
    % Fall contains only every tenth frame; use matching entries for the heatmap.
    state_id_Fall = state_id(mod(state_id, 10) == 0) / 10;
    state_id_real = state_id * step + start_frame - 1;

    %% 9. Compute house occupancy and stationary direction/axis preference
    CenXY_state = CenXY_plot(state_id, :);
    CenXY_in = zeros(length(CenXY_state), 8);  % stationary frames in houses.
    CenXY_in_all = zeros(length(CenXY_plot), 8); % all frames in houses.
    for j = 1:8
        CenXY_in(:, j) = inpolygon(CenXY_state(:, 1), CenXY_state(:, 2), Houses{j, 1}.Position(:, 1), resolution(1) - Houses{j, 1}.Position(:, 2));
        CenXY_in_all(:, j) = inpolygon(CenXY_plot(:, 1), CenXY_plot(:, 2), Houses{j, 1}.Position(:, 1), resolution(1) - Houses{j, 1}.Position(:, 2));
    end

    CenXY_out = 1 - sum(CenXY_in, 2); % stationary frames outside houses.

    % House order corresponds to 0, 45, ..., 315 degrees.
    Directions = 0:45:359;
    Counts = sum(CenXY_in);
    Norm_Counts = Counts / sum(Counts);
    [u_d, v_d] = pol2cart(deg2rad(Directions), Norm_Counts);
    img_dsi = u_d + sqrt(-1) * v_d;
    img_sum_dsi = sum(img_dsi);
    DSI = abs(img_sum_dsi);

    % Combine opposite houses and double their angles for axial preference.
    Orientations = [0 90 180 270];
    Counts = reshape(Counts, [4, 2]);
    Counts = sum(Counts, 2)';
    Norm_Counts = Counts / sum(Counts);
    [u_o, v_o] = pol2cart(deg2rad(Orientations), Norm_Counts);
    img_osi = u_o + sqrt(-1) * v_o;
    img_sum_osi = sum(img_osi);
    OSI = abs(img_sum_osi);

    RRE = length(state_id) / length(CenXY_plot);
    RIH = sum(sum(CenXY_in)) / length(CenXY_plot);

    R.u_d = u_d;
    R.v_d = v_d;
    R.img_dsi = img_dsi;
    R.DSI = DSI;
    R.u_o = u_o;
    R.v_o = v_o;
    R.img_osi = img_osi;
    R.OSI = OSI;
    R.RIH = RIH;
    R.RRE = RRE;
    R.dis = dis_trace;
    R.CenXY_in = CenXY_in;

    %% 10. Draw the summary figure
    H1 = figure('visible', 'off', 'position', [270, 95, 845, 845]);
    titlesize = 13;
    labelsize = 10;

    cl2 = rgb('Red'); % stationary related
    cl3 = rgb('Black'); % all frames
    clhouse = {'DodgerBlue', 'Aqua', 'LimeGreen', 'GreenYellow', 'Sandybrown', 'Red', 'Fuchsia', 'BlueViolet'};
    % House polygons.
    h1_1 = subplot(4, 3, 7);
    imagesc(F_averaged);
    hold on;
    colormap(h1_1, 'gray');
    axis off;
    for i = 1:8
        fill(Houses{i, 1}.Position(:, 1), Houses{i, 1}.Position(:, 2), rgb(clhouse{i}), 'FaceAlpha', .3, 'EdgeColor', 'none');
    end
    title({'Polygon of houses'; ' '}, 'fontsize', titlesize, 'color', 'k');
    hold off;
    scale = 100; % 0.7mm per pixel, 100 pixel/ 7cm
    % Displacement trace with stationary periods marked by house.
    subplot(4, 1, 1);
    plot(dis_trace, 'color', rgb('silver'), 'linewidth', 0.1);
    hold on;
    plot([0 length(CenXY_plot)], [threshold threshold], '-', 'color', cl2, 'linewidth', .1);
    ylim([-scale scale * 10.5]);
    box off;
    yy = zeros(length(state_id), 3);
    yy(:, 1) = -scale / 2;
    yy(:, 2) = scale;
    yy(:, 3) = NaN;

    for i = 1:8
        state_id_house = state_id(logical(CenXY_in(:, i)));
        h2_1 = bar(state_id_house, -yy(logical(CenXY_in(:, i)), 2), 'facecolor', rgb(clhouse{i}), 'edgecolor', 'none', 'showbaseline', 'off');
        alpha(h2_1, .4);
    end
    h2_1 = bar(state_id, -yy(:, 2), 'facecolor', cl3, 'edgecolor', 'none', 'showbaseline', 'off');
    alpha(h2_1, .02);

    xlabel('Time (hrs)');
    ylabel('Distance (cm)');
    legend({'Frames interval distance', 'Stationary threshold', 'B', 'BR', 'R', 'FR', 'F', 'FL', 'L', 'BL'}, 'box', 'off', 'fontsize', labelsize, 'NumColumns', 5);

    nstep = 3600 / (step * timeR); % how many datapoints for 1 hour

    set(gca, 'xlim', [0 length(CenXY_plot)], 'xtick', 0:nstep:length(CenXY_plot), 'xticklabel', 0:ceil(length(CenXY_plot) / nstep), ...
        'ytick', 0:scale:8 * scale, 'yticklabel', (0:8) * 7, 'fontsize', labelsize, 'Layer', 'top', 'TickLength', [0.006 0.015]);
    title({Mid; ['Video overview: minimum sleep section: 450secs, threshold: ', num2str(threshold * 70 / scale), 'mm']}, 'fontsize', titlesize, 'color', 'k');

    text(50, 8.9 * scale, ['Stationary(In-House): ', num2str(round(RIH, 2) * 100), '%'], 'fontsize', labelsize, 'HorizontalAlignment', 'left');
    text(50, 9.65 * scale, ['Stationary: ', num2str(round(RRE, 2) * 100), '%'], 'fontsize', labelsize, 'HorizontalAlignment', 'left');
    set(gca, 'tickdir', 'out');
    hold off;

    % Cumulative movement and occupancy of all eight houses.
    subplot(4, 1, 2);
    Counts_in_all = sum(CenXY_in_all);

    PI = Counts_in_all(1:8) / sum(Counts_in_all(1:8));
    In_house_ratio = sum(Counts_in_all(1:8)) / length(CenXY_plot);
    R.PI = PI;
    R.In_house_ratio = In_house_ratio;

    plot(cumsum(dis_trace), 'Color', rgb('steelblue'), 'LineWidth', 1);
    hold on;

    ylim([-2500 * scale 5000 * scale]);
    box off;
    yy = zeros(length(state_id), 3);
    yy(:, 1) = 0;
    yy(:, 2) = 4000 * scale;
    yy(:, 3) = NaN;
    yy2 = zeros(length(all_id), 1);
    yy2(:, 1) = -300 * scale;

    for i = 1:8
        plot(all_id(logical(CenXY_in_all(:, i))), i * yy2(logical(CenXY_in_all(:, i))), '|', 'color', [rgb(clhouse{i}), 0.4], 'markersize', 1);
    end

    for i = 1:8
        state_id_house = state_id(logical(CenXY_in(:, i)));
        h2_1 = bar(state_id_house, yy(logical(CenXY_in(:, i)), 2), 'facecolor', rgb(clhouse{i}), 'edgecolor', 'none', 'showbaseline', 'off');
        alpha(h2_1, .4);
    end
    h2_1 = bar(state_id, yy(:, 2), 'facecolor', cl3, 'edgecolor', 'none', 'showbaseline', 'off');
    alpha(h2_1, .02);

    xlabel('Time (hrs)');
    ylabel('Cumulative Distance (km)');
    set(gca, 'xlim', [0 length(CenXY_plot)], 'xtick', 0:nstep:length(CenXY_plot), 'xticklabel', 0:ceil(length(CenXY_plot) / nstep), ...
        'ytick', (0:100000:400000) * scale / 100, 'yticklabel', (0:4) * .07, 'fontsize', labelsize, 'Layer', 'top', 'TickLength', [0.006 0.015]);
    text(50, 4500 * scale, ['Minimum distance traveled overall: ', num2str(round((7 / scale) * sum(dis_trace) / 100000, 2)), ' km'], 'fontsize', labelsize, 'HorizontalAlignment', 'left');
    text(50, 4000 * scale, ['In house ratio: ', num2str(round(In_house_ratio, 2))], 'fontsize', labelsize, 'HorizontalAlignment', 'left');
    text(50, 3500 * scale, ['Preference index: ', strjoin(compose('%.2f', round(PI, 2)))], 'fontsize', labelsize, 'HorizontalAlignment', 'left');
    set(gca, 'tickdir', 'out');
    title({'Distance traveled & house visiting plot'}, 'fontsize', titlesize, 'color', 'k');
    hold off;
    % Heatmaps use the stored masks, not every analyzed frame.
    h5_1 = subplot(4, 3, 8);
    imagesc(squeeze(sum(Fall, 1)));
    xlim([0 resolution(2)]);
    ylim([0 resolution(1)]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    title({'Heatmap of all frames'; ' '}, 'fontsize', titlesize, 'color', 'k');
    colormap(h5_1, 'turbo');

    h7_1 = subplot(4, 3, 9);
    imagesc(squeeze(sum(Fall(state_id_Fall, :, :), 1)));
    xlim([0 resolution(2)]);
    ylim([0 resolution(1)]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    title({'Heatmap of stationary frames'; ' '}, 'fontsize', titlesize);
    colormap(h7_1, 'turbo');

    % Spatial positions and stationary counts for each house.
    subplot(4, 3, 10);
    plot(CenXY_plot(:, 1), CenXY_plot(:, 2), '.', 'markersize', 1, 'color', rgb('gray'));
    hold on;
    xlim([0 resolution(2)]);
    ylim([0 resolution(1)]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    Counts = sum(CenXY_in);

    for i = 1:8
        plot(CenXY_state(logical(CenXY_in(:, i)), 1), CenXY_state(logical(CenXY_in(:, i)), 2), '.', 'color', rgb(clhouse{i}), 'markersize', 2);
        xlim([0 resolution(2)]);
        ylim([0 resolution(1)]);
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
    plot(CenXY_state(logical(CenXY_out), 1), CenXY_state(logical(CenXY_out), 2), '+', 'color', cl3, 'markersize', 3);

    xC = [580 960 1150 960 560 270 60 270] * resolution(2) / 1280;
    yC = [1000 850 500 150 24 150 500 850] * resolution(1) / 1024;

    for i = 1:8
        text(xC(i), yC(i), num2str(Counts(i)), 'fontsize', labelsize);
    end
    title({'In-house sleeping frames'; ' '}, 'fontsize', titlesize);
    box off;
    axis off;
    hold off;

    % Directional vector plot.
    subplot(4, 3, 11);
    pd = angle(img_sum_dsi);

    pd_8 = deg2rad(Directions([1:end, 1]));
    ln_8 = abs([img_dsi, img_dsi(1)]);
    for i = 1:50
        if i < 50
            polarplot(pd_8, ln_8 * i / 50, 'color', [cl3 .1], 'linewidth', 0.05);
            hold on;
        else
            polarplot(pd_8, ln_8 * i / 50, 'color', [cl3 .5], 'linewidth', 0.5);
            hold on;
        end
    end
    polarplot([0 pd], [0, DSI], 'color', [cl3 .9], 'linewidth', 1.5);

    axp = gca;
    axp.Layer = 'bottom';
    axp.ThetaDir = 'clockwise';
    axp.ThetaZeroLocation = 'Top';
    axp.ThetaGrid = 'on';
    axp.RTick = [0 0.5 1];
    axp.RLim = [0 1];
    axp.FontSize = labelsize;
    axp.ThetaTick = [0 90 180 270];
    axp.ThetaTickLabel = {'Back'; 'Right'; 'Front'; 'Left'};
    axp.LineWidth = 0.5;
    axp.GridAlpha = 0.15;
    title(['Direction preference: ', num2str(round(DSI, 2))], 'fontsize', titlesize);
    hold off;
    R.pd = pd;
    % Axial vector plot.
    subplot(4, 3, 12);
    po = angle(sum(img_osi));
    po_4 = deg2rad([0 90 180 270 0]);
    ln_4 = abs([img_osi, img_osi(1)]);

    for i = 1:50
        if i < 50
            polarplot(po_4, ln_4 * i / 50, 'color', [cl3 .1], 'linewidth', 0.05);
            hold on;
        else
            polarplot(po_4, ln_4 * i / 50, 'color', [cl3 .5], 'linewidth', 0.5);
            hold on;
        end
    end
    polarplot([0 po], [0, OSI], 'color', [cl3 .9], 'linewidth', 1.5);

    axp = gca;
    axp.Layer = 'bottom';
    axp.ThetaDir = 'clockwise';
    axp.ThetaZeroLocation = 'Top';
    axp.ThetaTick = [0 90 180 270];
    axp.ThetaTickLabel = {'B-F'; 'BR-FL'; 'R-L'; 'FR-BL'};
    axp.ThetaGrid = 'on';
    axp.RTick = [0 .5 1];
    axp.RLim = [0 1];
    axp.FontSize = 10;
    axp.LineWidth = 0.5;
    axp.GridAlpha = 0.15;

    R.po = po;
    title(['Axis preference: ', num2str(round(OSI, 2))], 'fontsize', titlesize);
    hold off;

    %% 11. Save summary figures and analysis variables
    orient landscape;
    set(H1, 'papersize', [9 12], 'paperposition', [.01 .01 8.98 11.98]);

    fstd_str = '450secs';
    % Preserve the original output names; this label describes the default threshold.
    threshold_str = '7mm';

    if full == 0
        saveas(H1, ['Summary_', Mid0, '_', fstd_str, '_', threshold_str, '_full.pdf']);
        saveas(H1, ['Summary_', Mid0, '_', fstd_str, '_', threshold_str, '_full.png']);
    else
        saveas(H1, ['Summary_', Mid0, '_', fstd_str, '_', threshold_str, '.pdf']);
        saveas(H1, ['Summary_', Mid0, '_', fstd_str, '_', threshold_str, '.png']);
    end

    if full == 0
        save(fullfile(path2, ['Results_', fstd_str, '_', threshold_str, '_full.mat']), 'CenXY_plot', 'CenXY_state', 'CenXY_in', 'state_id_Fall', 'Fall', 'R', 'Houses', 'state_id_real');
    else
        save(fullfile(path2, ['Results_', fstd_str, '_', threshold_str, '.mat']), 'CenXY_plot', 'CenXY_state', 'CenXY_in', 'state_id_Fall', 'Fall', 'R', 'Houses', 'state_id_real');
    end

    close all;
end
