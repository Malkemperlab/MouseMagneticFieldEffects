function R = HouseChoosingTest(data, tl, VecID, MouseID, ColorID, DircomID)
    % HOUSECHOOSINGTEST Plot house-choice directions and circular statistics.
    %   R = HouseChoosingTest(data, tl, VecID, MouseID, ColorID, DircomID)
    %
    % Required inputs
    %   data      Column vector of already-shifted individual angles, in degrees.
    %             The desired reference direction must correspond to 0 degrees.
    %   tl        Title identifier used as the output filename prefix.
    %
    % Optional inputs (defaults in brackets)
    %   VecID     0: house-choice marker layout and grouping correction [0].
    %             1: vector-summary marker layout without the correction.
    %             Both modes retain the original rounding to 45-degree bins.
    %   MouseID   Optional labels in the original data order [empty].
    %   ColorID   Marker color indices, 1 to 8 [ones(length(data),1)].
    %             Marker edges use the same color, except color 8 uses color 1.
    %   DircomID  1: label the zero-degree reference tNorth [1].
    %             2: label the zero-degree reference mHome (trained direction).
    %             This option changes the label only; tN is fixed at 0.
    %
    % Output
    %   R = [r_direction, p_direction]. The r value includes the original
    %   1.0262 grouping correction when VecID is 0. The plotted vector length
    %   uses the uncorrected value. The directional Rayleigh p-value formula
    %   and confidence-interval drawing are retained.
    %
    % Files: <tl>_DP.svg and <tl>_DP.png in the current folder.
    % Requires CircStat functions circ_r and circ_mean, plus rgb.m.
    % The function creates hidden figures and closes figures after saving.
    % Author: Li Zhang. Original version: 2022-03-08.

    %% 1. Plot palette and optional arguments
    tN = 0; % fixed reference direction for the already-shifted data
    % define color order of datapoint sets
    color_dp = {'slategray', 'DodgerBlue', 'MediumSeaGreen', 'salmon', 'LightSkyBlue', 'Orchid', 'orange', 'silver'};
    % Defaults for the optional plotting arguments.
    if nargin < 3
        VecID = 0;
    end
    if nargin < 4
        MouseID = [];
    end
    if nargin < 5
        ColorID = ones(length(data), 1);
    end
    if nargin < 6
        DircomID = 1;
    end
    if ~isscalar(DircomID) || ~ismember(DircomID, [1, 2])
        error('HouseChoosingTest:InvalidDircomID', ...
            'DircomID must be 1 (tNorth) or 2 (mHome).');
    end

    %% 2. Wrap the supplied angles and assign 45-degree bins
    disp(tl);
    disp(['number of datapoints: ', num2str(length(data))]);
    % The caller has aligned the data; wrap into [0, 360) without rotating it.
    data_shift = mod(data, 360);

    % round data to nearest 45 degrees for binning
    roundN = 45; % default 45
    po = [0:roundN:360];
    po = po(2:end);
    order = [1:length(po)];
    data_shift = round(data_shift / roundN) * roundN;
    data_shift(data_shift == 0) = 360;

    %% 3. Directional preference: circular mean and Rayleigh test
    figure('position', [50 50 666 500], 'visible', 'off');
    ln = [];
    % count datapoints per angular bin
    for i = 1:length(order)
        ln(order(i), 1) = length(find(data_shift == po(i)));
    end
    M = find(ln > 0);
    [~, I] = sort(data_shift);

    % compute mean resultant length and mean angle using circ toolbox
    r = circ_r(deg2rad(data_shift(:)));
    [sita, ul, ll] = circ_mean(deg2rad(data_shift(:))); % sita: mean angle;
    r1 = r;

    % correction factor for grouping
    if VecID == 0
        r1 = r1 * 1.0262; % group = 8, factor = 1.0262;
    end

    % compute Rayleigh's test p-value
    n = length(data_shift);
    rayleigh_R = n * r1;
    pval = exp(sqrt(1 + 4 * n + 4 * (n^2 - rayleigh_R^2)) - (1 + 2 * n));
    p1 = pval;

    % pre-draw empty polar rings for layout
    for i = M
        for j = 1:5
            polarscatter(deg2rad(po(i)), 1.4 + j * .28, 85, 'w');
            hold on;
        end
    end

    % draw resultant vector arrow
    polarplot([0 sita], [0 r * 1.5], 'k', 'linewidth', 1);
    polarplot([sita + deg2rad(8) sita], [r * 1.5 - .1 r * 1.5], 'k', 'linewidth', 1);
    polarplot([sita - deg2rad(8) sita], [r * 1.5 - .1 r * 1.5], 'k', 'linewidth', 1);

    % if significant, plot confidence interval
    if p1 < 0.05
        polarplot([0 ul 0 ll], [0 1.4 0 1.4], '--k', 'linewidth', 1);
    end

    % adjust polar plot appearance
    if VecID == 1
        set(gca, 'RTick', [0 1.5], 'RTickLabel', [], 'RLim', [0 4], 'ThetaZeroLocation', 'Top', ...
            'ThetaDir', 'clockwise', 'fontsize', 18);
    else
        set(gca, 'RTick', [0 1.5], 'RTickLabel', [], 'RLim', [0 3], 'ThetaZeroLocation', 'Top', ...
            'ThetaDir', 'clockwise', 'fontsize', 18);
    end

    % scatter datapoints in rings based on bin frequency
    dpo = 8; % degree shift for spacing
    mid = 0;
    for i = 1:length(M)
        if ln(M(i), 1) > 0
            if VecID == 0
                n = ln(M(i), 1);
                ks = mod(n, 3);
                kn = floor(n / 3);
                % generate pattern for multi-ring scatter arrangement
                if ks == 2
                    poks = [po(M(i)) - dpo, po(M(i)), po(M(i)) + dpo];
                    poks = [repmat(poks, 1, kn), po(M(i)) - dpo / 2, po(M(i)) + dpo / 2];
                elseif ks == 1
                    poks = [po(M(i)) - dpo, po(M(i)), po(M(i)) + dpo];
                    poks = [repmat(poks, 1, kn), po(M(i))];
                elseif ks == 0
                    poks = [po(M(i)) - dpo, po(M(i)), po(M(i)) + dpo];
                    poks = [repmat(poks, 1, kn)];
                end
                % plot individual datapoints
                for j = 1:n
                    mid = mid + 1;
                    if ColorID(I(mid)) == 8
                        polarscatter(deg2rad(poks(j)), 1.4 + ceil(j / 3) * .22, 32 + ceil(j / 3) * 7, rgb(color_dp{ColorID(I(mid))}), 'filled', 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', rgb(color_dp{1}));
                    else
                        polarscatter(deg2rad(poks(j)), 1.4 + ceil(j / 3) * .22, 32 + ceil(j / 3) * 7, rgb(color_dp{ColorID(I(mid))}), 'filled', 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', rgb(color_dp{ColorID(I(mid))}));
                    end
                    % optionally annotate MouseID
                    if ~isempty(MouseID)
                        text(deg2rad(poks(j)), 1.4 + ceil(j / 3) * .22, MouseID(I(mid)), 'fontsize', 3, 'HorizontalAlignment', 'center', 'color', 'white');
                    end
                end
            elseif VecID == 1
                n = ln(M(i), 1);
                for j = 1:n
                    mid = mid + 1;
                    if ColorID(I(mid)) == 8
                        polarscatter(deg2rad(po(M(i))), 1.4 + j * .24, 35, rgb(color_dp{ColorID(I(mid))}), 'filled', 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', rgb(color_dp{1}));
                    else
                        polarscatter(deg2rad(po(M(i))), 1.4 + j * .24, 35, rgb(color_dp{ColorID(I(mid))}), 'filled', 'MarkerFaceAlpha', .5, 'MarkerEdgeColor', rgb(color_dp{ColorID(I(mid))}));
                    end
                    if ~isempty(MouseID)
                        text(deg2rad(po(M(i))), 1.4 + j * .24, MouseID(I(mid)), 'fontsize', 3, 'HorizontalAlignment', 'center', 'color', 'white');
                    end
                end
            end
        end
    end

    % label direction ticks depending on DircomID
    if DircomID == 1
        set(gca, 'ThetaTick', tN, 'ThetaTickLabel', 'tNorth', 'RColor', [0 0 0], ...
            'GridAlpha', 1, 'ThetaGrid', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0]);
    elseif DircomID == 2
        set(gca, 'ThetaTick', tN, 'ThetaTickLabel', 'mHome', 'RColor', [0 0 0], ...
            'GridAlpha', 1, 'ThetaGrid', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0]);
    end

    % add title with results
    title({'Direction Preference'; ['r=', num2str(round(r1, 2)), ' p=', num2str(round(pval, 7))]}, 'fontsize', 15, 'fontweight', 'normal');
    hold off;
    % Save the directional plot in editable and raster formats.
    orient landscape;
    set(gcf, 'papersize', [3.5 3.5], 'paperposition', [0 0 3.5 3.5]);
    saveas(gcf, [tl, '_DP.svg']);
    saveas(gcf, [tl, '_DP.png']);
    close all;

    %% 4. Return directional statistics
    R = [r1, p1];
end
