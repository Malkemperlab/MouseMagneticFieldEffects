%% Run_HouseChoosingTest
% Generate two polar plots for every Panel_id in circular_datasets.xlsx:
%   Left:  trained magnetic home direction at zero, labelled mHome.
%   Right: geographic north at zero, labelled tNorth.
%
% SETUP
% 1. Save this script and circular_datasets.xlsx in the same folder.
% 2. Put HouseChoosingTest.m in that folder, or change functionFolder 
%    below to its location.
% 3. Add CircStat (circ_r and circ_mean) and rgb.m to the MATLAB path.
% 4. Run this whole script. Adjust the settings in section 1 if needed.
%
% HOW THE TRAINED-MAGNETIC-NORTH COLUMN WAS CALCULATED
% "First_Visit_House (aligned to TrainedMagNorth)" was calculated as:
%
%   first visit house aligned to geoNorth
%       - testing mN
%       - trained direction deviation from mN
%
% In variable form, with all angles in degrees:
%   angle_mHome = angle_geoNorth - testing_mN - trained_deviation_from_mN;
%
% testing_mN is the magnetic north direction during testing, expressed
% relative to geographic north. trained_deviation_from_mN is the trained
% home direction's signed angular offset from magnetic north.
% Use the same angular sign convention for all three terms.
%
% Excel column B already contains this calculated angle; column C contains
% the geographic-north angle. This script reads both columns directly.
% It does not subtract testing_mN or the trained offset again. These source
% offsets are not separate columns in the supplied workbook.
% Negative angles and angles outside 0..360 are valid: HouseChoosingTest
% applies mod(data,360), then retains its rounding to 45-degree bins.
% DircomID only selects the displayed label; it does not align the angles.
%
% OUTPUTS (inside outputFolder)
%   <panel>_mHome_DP.svg / .png   Individual mHome polar plot.
%   <panel>_tNorth_DP.svg / .png  Individual tNorth polar plot.
%   <panel>_paired.png           Both plots side by side, mHome on the left.
%   HouseChoosingTest_statistics.csv   N, resultant length r, and p for both.
%
% The original function saves and closes its own figures. Therefore, the
% paired PNG is assembled from its two exported PNGs after both calls finish.
% Use the individual SVG files when you need editable vector graphics.

%% 1. Paths and plotting settings
scriptFolder = 'xxx:\HouseChoosingTest_analysis'; % replace with the path contains data and scripts
if isempty(scriptFolder)
    scriptFolder = pwd;
end

excelFile = fullfile(scriptFolder, 'circular_datasets.xlsx');
sheetName = 'Sheet1';
functionFolder = scriptFolder;
outputFolder = fullfile(scriptFolder, 'HouseChoosingTest_plots');

VecID = 0;          % House-choice layout, including the function's grouping correction.
defaultColorID = 1; % Marker color: 1..8, using the palette in HouseChoosingTest.
% No mouse IDs are supplied because the workbook has no MouseID column.

% If needed, uncomment and edit these paths:
% addpath('C:\MATLAB\CircStat');
% addpath('C:\MATLAB\ColorFunctions'); % Folder containing rgb.m.

%% 2. Make the supplied plotting function available
functionFile = fullfile(functionFolder, 'HouseChoosingTest.m');
if ~isfile(functionFile)
    error('RunHouseChoosingTest:MissingFunction', ...
        ['Save the supplied function as HouseChoosingTest.m in functionFolder. ' ...
         'Its downloaded name may be HouseChoosingTest.m.']);
end
addpath(functionFolder, '-begin');

requiredFunctions = {'HouseChoosingTest', 'circ_r', 'circ_mean', 'rgb'};
for k = 1:numel(requiredFunctions)
    if exist(requiredFunctions{k}, 'file') == 0
        error('RunHouseChoosingTest:MissingDependency', ...
            'Add the folder containing %s.m to the MATLAB path.', requiredFunctions{k});
    end
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% 3. Read panel labels and the two already-aligned angle columns
T = readtable(excelFile, 'Sheet', sheetName, ...
    'VariableNamingRule', 'preserve', 'TextType', 'string');

% Normalize repeated spaces or line breaks in the Excel headers only.
headers = strtrim(regexprep(string(T.Properties.VariableNames), '\s+', ' '));
expectedHeaders = ["Panel_id", ...
    "First_Visit_House (aligned to TrainedMagNorth)", ...
    "First_Visit_House (aligned to geoNorth)"];
[found, columnIndex] = ismember(expectedHeaders, headers);
if ~all(found)
    error('RunHouseChoosingTest:MissingColumn', ...
        'Missing Excel column(s): %s', char(strjoin(expectedHeaders(~found), ', ')));
end

panelIDs = strtrim(string(T{:, columnIndex(1)}));
angles_mHome = str2double(string(T{:, columnIndex(2)}));
angles_tNorth = str2double(string(T{:, columnIndex(3)}));

% Require both angles for each row so the two plots use the same animals.
% Zero is a valid angle. Missing or nonnumeric values stop the script instead
% of silently excluding observations from either alignment.
invalidRows = ismissing(panelIDs) | strlength(panelIDs) == 0 ...
    | isnan(angles_mHome) | isinf(angles_mHome) ...
    | isnan(angles_tNorth) | isinf(angles_tNorth);
if any(invalidRows)
    excelRows = find(invalidRows) + 1; % Add one for the Excel header row.
    error('RunHouseChoosingTest:InvalidData', ...
        'Each row needs Panel_id and two finite angles. Check Excel row(s): %s', ...
        sprintf('%d ', excelRows));
end
if isempty(panelIDs)
    error('RunHouseChoosingTest:EmptyTable', 'The selected worksheet has no data rows.');
end

% Process every panel automatically, keeping its order in the Excel table.
panels = unique(panelIDs, 'stable');
nPanels = numel(panels);
N = zeros(nPanels, 1);
R_mHome = nan(nPanels, 2);
R_tNorth = nan(nPanels, 2);

%% 4. Generate the two polar plots and the paired image for each panel
for k = 1:nPanels
    thisPanel = panels(k);
    rows = panelIDs == thisPanel;
    N(k) = sum(rows);
    colors = repmat(defaultColorID, N(k), 1);

    % Use the same row selection and marker colors for both alignments.
    data_mHome = angles_mHome(rows);
    data_tNorth = angles_tNorth(rows);

    % Example filename conversion: Fig.4D becomes Fig_4D.
    panelTag = regexprep(char(thisPanel), '[^A-Za-z0-9_-]', '_');
    prefix_mHome = fullfile(outputFolder, [panelTag, '_mHome']);
    prefix_tNorth = fullfile(outputFolder, [panelTag, '_tNorth']);
    fprintf('\n%s: %d observations, generating mHome and tNorth plots.\n', ...
        char(thisPanel), N(k));

    % The final argument controls only the reference label:
    %   2 = mHome, paired with the TrainedMagNorth column.
    %   1 = tNorth, paired with the geoNorth column.
    R_mHome(k, :) = HouseChoosingTest(data_mHome, prefix_mHome, ...
        VecID, [], colors, 2);
    R_tNorth(k, :) = HouseChoosingTest(data_tNorth, prefix_tNorth, ...
        VecID, [], colors, 1);

    % Create the paired figure after both function calls, because each call
    % closes figures. The original individual SVG and PNG files are retained.
    pairedFigure = figure('Visible', 'off', 'Color', 'white', ...
        'Position', [100, 100, 1200, 600]);
    subplot(1, 2, 1);
    image(imread([prefix_mHome, '_DP.png']));
    axis image off;
    title('Trained magnetic direction: mHome', 'FontSize', 13);

    subplot(1, 2, 2);
    image(imread([prefix_tNorth, '_DP.png']));
    axis image off;
    title('Geographic north: tNorth', 'FontSize', 13);
    sgtitle(char(thisPanel), 'Interpreter', 'none', 'FontSize', 16);

    set(pairedFigure, 'PaperPositionMode', 'auto');
    print(pairedFigure, fullfile(outputFolder, [panelTag, '_paired.png']), ...
        '-dpng', '-r300');
    close(pairedFigure);
end

%% 5. Save the statistics returned by HouseChoosingTest
% The function returns [r, p]. Keep its calculations and correction unchanged.
results = table(panels, N, R_mHome(:, 1), R_mHome(:, 2), ...
    R_tNorth(:, 1), R_tNorth(:, 2), 'VariableNames', ...
    {'Panel_id', 'N', 'r_mHome', 'p_mHome', 'r_tNorth', 'p_tNorth'});
writetable(results, fullfile(outputFolder, 'HouseChoosingTest_statistics.csv'));
disp(results);
fprintf('\nFinished. Plots and statistics were saved in:\n%s\n', outputFolder);

