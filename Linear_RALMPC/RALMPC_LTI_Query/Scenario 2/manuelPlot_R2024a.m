clear
clc
close all
%search for all datasets
mat = dir('/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Data/*.mat'); 
%% Loop through all data sets
Numberofiteration=[];
for q = 1:length(mat) 
    load(mat(q).folder+"/"+mat(q).name);
    X_RALMPC_N{q}=X_RALMPC;
    J_RALMPC_N_all{q}=J_RALMPC;
    J_RALMPC_N{q}=J_RALMPC{end};
    Numberofiteration(q)=N;
    T_cpu_RALMPC{q}=mean(vertcat(t_cpu_RALMPC{:}));
end
T_cpu_RAMPC=mean(vertcat(t_cpu_RAMPC{:}));
[id id] = sort(Numberofiteration,"descend");
Numberofiteration=Numberofiteration(id);
k=1;
for i=id
X_RALMPC_N_{k}=X_RALMPC_N{i};
J_RALMPC_N_{k}=J_RALMPC_N{i};
T_cpu_RALMPC_{k}=T_cpu_RALMPC{i};
J_RALMPC_N_all_{k}=J_RALMPC_N_all{i};
k=k+1;
end
X_RALMPC_N=X_RALMPC_N_;
J_RALMPC_N=J_RALMPC_N_;
T_cpu_RALMPC=T_cpu_RALMPC_;
J_RALMPC_N_all=J_RALMPC_N_all_;
%% Plot 1: Comparison Computation time over N
% MATLAB script for a bar plot with colorblind-friendly colors and legend

% Ensure valid indices and values
[isMember, indices] = ismember([18,12,8,6], Numberofiteration);

% Data for the plot
values = [T_cpu_RAMPC, cell2mat(T_cpu_RALMPC(indices))]; % Combine reference and comparison values

% Create corresponding category labels for the x-axis
categories = cellstr(['N=25', compose("N=%d", Numberofiteration(indices))]);

% Colorblind-friendly colors
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8)
compColor = [230, 159, 0] / 255;  % Orange (#E69F00)

% Create the bar plot
figure;
b = bar(values, 'FaceColor', 'flat'); % Plot bars without categorical alignment

% Assign colors: First bar (RAMPC) different, others (RALMPC) the same
b.CData(1, :) = refColor; % RAMPC bar (reference)
b.CData(2:end, :) = repmat(compColor, length(values)-1, 1); % RALMPC bars (comparison)

% Manually set x-axis labels
xticks(1:length(categories)); % Set the ticks based on number of categories
xticklabels(categories); % Assign the dynamic labels to the ticks

% Labels and formatting
ylabel('Computation Time [s]', 'Interpreter', 'latex', 'FontSize', 14);

grid on;

% Axis settings
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Create the legend manually with colored patches
hold on;
h1 = patch(nan, nan, refColor); % Create a patch for RAMPC
h2 = patch(nan, nan, compColor); % Create a patch for RALMPC
legend([h1, h2], {'RAMPC', 'RALMPC'}, 'Location', 'best', 'Interpreter', 'latex');

% Export as PDF and EPS for LaTeX
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_ComputationTime.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_ComputationTime.eps');

disp('PaperPlot_ComputationTime.pdf and PaperPlot_ComputationTime.eps');
%% Plot 2: Comparison Cost time over N
% MATLAB script for a bar plot with colorblind-friendly colors and legend

% Ensure valid indices and values
[isMember, indices] = ismember([18,12,8,6], Numberofiteration);

% Data for the plot
values = [cell2mat(J_RAMPC(end)), cell2mat(J_RALMPC_N(indices))]; % Combine reference and comparison values

% Create corresponding category labels for the x-axis
categories = cellstr(['N=25', compose("N=%d", Numberofiteration(indices))]);

% Colorblind-friendly colors
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8)
compColor = [230, 159, 0] / 255;  % Orange (#E69F00)

% Create the bar plot
figure;
b = bar(values, 'FaceColor', 'flat'); % Plot bars without categorical alignment

% Assign colors: First bar (RAMPC) different, others (RALMPC) the same
b.CData(1, :) = refColor; % RAMPC bar (reference)
b.CData(2:end, :) = repmat(compColor, length(values)-1, 1); % RALMPC bars (comparison)

% Manually set x-axis labels
xticks(1:length(categories)); % Set the ticks based on number of categories
xticklabels(categories); % Assign the dynamic labels to the ticks

% Add a horizontal line for the optimal cost (colorblind-friendly)
yline(J_OS, '--', 'Optimal Cost', 'Color', '#0072B2', ...
    'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);

% Labels and formatting
ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14);
ylim([100,160]);

grid on;

% Axis settings
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Create the legend manually with colored patches
hold on;
h1 = patch(nan, nan, refColor); % Create a patch for RAMPC
h2 = patch(nan, nan, compColor); % Create a patch for RALMPC
legend([h1, h2], {'RAMPC', 'RALMPC'}, 'Location', 'northeast', 'Interpreter', 'latex');

% Export as PDF and EPS for LaTeX
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Cost.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Cost.eps');

disp('PaperPlot_ComputationTime.pdf and PaperPlot_ComputationTime.eps');
%% Plot 3: x1 over time
% MATLAB script with improved colorblind-friendly colors and legend handling

% Find the index of iteration N=8
index = find(Numberofiteration == 8);

% Extract the state for RALMPC and RAMPC iterations
X_RALMPC_iteration = X_RALMPC_N{index};

% Get the last state of RAMPC iteration
state_RAMPC_iteration = X_RAMPC{end};

% Time vector
time = 0:T_s:T_s*(length(state_RAMPC_iteration)-1);

% Define colorblind-friendly colors (Color Universal Design - CUD)
X_Init_Color = [0.2, 0.2, 0.8];    % Dark blue for initial trajectory
RAMPC_Color = [0.0, 0.45, 0.7];  % Blue (RAMPC)
OS_Color = [0.8, 0.4, 0.0];      % Orange (OS)
RALMPC_First_Color = [0.6, 0.8, 0.2]; % Light Green (RALMPC First Iteration)
RALMPC_Last_Color = [0.8, 0.17, 0.55]; % Magenta (RALMPC Last Iteration)

% Create the figure
figure(3);
hold on;

% Plot initial trajectory (X_bar_initial)
plot(time, X_bar_inital(1,:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'Color', X_Init_Color, 'DisplayName', 'Initial Trajectory');

% Plot RAMPC (solid blue line)
plot(time, state_RAMPC_iteration(1,:), '--', 'LineWidth', 2, 'DisplayName', 'RAMPC', 'Color', RAMPC_Color);

% Plot OS (solid line with markers)
plot(time, X_OS(1,:), '-*', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'OS', 'Color', OS_Color);

% Plot RALMPC (first iteration - light green, NO legend entry)
state_RALMPC_iteration = X_RALMPC_iteration{1};
plot(time, state_RALMPC_iteration(1,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (First It.)','Color', RALMPC_First_Color); % No 'DisplayName'

% Plot RALMPC (last iteration - magenta, legend included)
state_RALMPC_iteration = X_RALMPC_iteration{end};
plot(time, state_RALMPC_iteration(1,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (Last It.)', 'Color', RALMPC_Last_Color);

% Add a horizontal line at y = 0 (EXCLUDED from legend)
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Add labels, title, and grid
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
xlim([0,5.9]);

% Set axis properties for better readability
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Add a legend 
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

% Hold off to stop adding more plots
hold off;
% Export as PDF and EPS for LaTeX
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_timex1.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_timex1.eps');

disp('PaperPlot_ComputationTime.pdf and PaperPlot_ComputationTime.eps');

%% Plot 4:
figure(4);
J_RALMPC_iteration = J_RALMPC_N_all{index}; % Extract cost values
J_values = cell2mat(J_RALMPC_iteration); % Convert to numerical array

% Define initial cost (assuming it is known or the first value)
initialCost = J_values(1); 

% Define optimal cost (assumed known, replace with actual value)
optimalCost = 140; % Replace with actual optimal cost

% Total number of bars
numBars = length(J_values);

% Select the first 4 and the last value
selectedIndices = [1:4, numBars]; % Take first 4 and last
selectedValues = J_values(selectedIndices);

% Add initial cost at x = 0
xPositions = [0, 1:4, 5.5, 7];  % Position 5.5 is where we place "..."
xLabels = ["Initial", string(1:4), "...", string(numBars)]; % Labels with dots in the gap
barValues = [J_0, selectedValues]; % Include initial cost

% Create a bar plot (excluding dots)
b = bar(xPositions([1, 2:5, 7]), barValues, 'FaceColor', 'flat'); % Enable color assignment
b.CData = compColor; % Assign colors to bars

% Adjust x-ticks and labels
xticks(xPositions);
xticklabels(xLabels);

% Set axis limits and grid
ylim([130, 155]);
grid on;

% Add a horizontal line for the optimal cost (colorblind-friendly)
yline(J_OS, '--', 'Optimal Cost', 'Color', '#0072B2', ...
    'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);

% Add labels and title
xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14);

% Improve appearance
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Export as PDF and EPS for LaTeX
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Costiteration.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Costiteration.eps');

disp('PaperPlot_Costiteration.pdf and PaperPlot_Costiteration.eps');

%% Figure 5: Improved 2D State Trajectory Plot for CDC
figure(5); clf; hold on;

% Find the index of iteration N=8
index = find(Numberofiteration == 8);

% Extract the state for RALMPC and RAMPC iterations
X_RALMPC_iteration = X_RALMPC_N{index};

% Get the last state of RAMPC iteration
state_RAMPC_iteration = X_RAMPC{end};

% Define colorblind-friendly colors
Polyhedron_Color = [0.8, 0.8, 1];  % Light blue for polyhedron (excluded from legend)
X_Init_Color = [0.2, 0.2, 0.8];    % Dark blue for initial trajectory

RAMPC_Color = [0.0, 0.45, 0.7];  % Blue (RAMPC)
OS_Color = [0.8, 0.4, 0.0];      % Orange (OS)
RALMPC_First_Color = [0.6, 0.8, 0.2]; % Light Green (RALMPC First Iteration)
RALMPC_Last_Color = [0.8, 0.17, 0.55]; % Magenta (RALMPC Last Iteration)

% Plot the Polyhedron constraints (EXCLUDED from the legend)
for k = 1:length(X_bar_inital)
    X_poly = Polyhedron(H, S_inital(k) * ones(size(H,1),1) + H * X_bar_inital(:,k));
    h = plot(X_poly, 'Color', Polyhedron_Color, 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
end

% Plot initial trajectory (X_bar_initial)
plot(X_bar_inital(1,:), X_bar_inital(2,:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'Color', X_Init_Color, 'DisplayName', 'Initial Trajectory');

% Plot RAMPC trajectory
plot(state_RAMPC_iteration(1,:), state_RAMPC_iteration(2,:), '--', 'LineWidth', 2, 'DisplayName', 'RAMPC', 'Color', RAMPC_Color);
 
% Plot OS trajectory
plot(X_OS(1,:), X_OS(2,:), '-*', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'OS', 'Color', OS_Color);

% Plot RALMPC first iteration (lighter color, NO legend entry)
plot(X_RALMPC_iteration{1}(1,:), X_RALMPC_iteration{1}(2,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (First It.)','Color', RALMPC_First_Color); % No 'DisplayName'
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude first iteration from legend

% Plot RALMPC last iteration (included in legend)
plot(X_RALMPC_iteration{end}(1,:), X_RALMPC_iteration{end}(2,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (Last It.)', 'Color', RALMPC_Last_Color);

% Define rectangle limits
x1_min = -0.2; x1_max = 4.1;
x2_min = -5; x2_max = 5;

% Compute width and height
width = x1_max - x1_min;
height = x2_max - x2_min;

% Draw the rectangle (black color, thicker line, not in legend)
rectangle('Position', [x1_min, x2_min, width, height], 'EdgeColor', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off');

% Set axis limits
xlim([-0.4, 4.3]); 
ylim([-5.3, 5.3]);
% Formatting
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
title('State Trajectory Comparison', 'Interpreter', 'latex', 'FontSize', 16);
grid on; 
%axis equal;

% Improve appearance
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Add a legend (Polyhedra and first RALMPC iteration excluded)
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

hold off;

% Export as PDF and EPS for LaTeX
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Costiteration.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Costiteration.eps');

disp('PaperPlot_ComputationTime.pdf and PaperPlot_ComputationTime.eps');

%% Combined Bar Plot: Computation Time & Cost (Side-by-Side, Black Y-Labels & Cost Limit)
figure;
hold on;

% Ensure valid indices and values
[isMember, indices] = ismember([18,12,8,6], Numberofiteration);

% Data for computation time (left y-axis)
compTimeValues = [T_cpu_RAMPC, cell2mat(T_cpu_RALMPC(indices))];

% Data for cost (right y-axis)
costValues = [cell2mat(J_RAMPC(end)), cell2mat(J_RALMPC_N(indices))];

% Create corresponding category labels for the x-axis
categories = cellstr(['N=25', compose("N=%d", Numberofiteration(indices))]);

% Colorblind-friendly colors
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8) for RAMPC
compColor = [230, 159, 0] / 255;  % Orange (#E69F00) for RALMPC

% Number of categories
numGroups = length(categories);
xPos = 1:numGroups; % X-axis positions

% ---- Plot Computation Time on Left Y-Axis ----
yyaxis left;
b1 = bar(xPos - 0.2, compTimeValues, 0.4, 'FaceColor', 'flat', 'EdgeColor', 'none'); % Left-shift bars
b1.CData(1, :) = refColor; % RAMPC (Computation Time) in blue
b1.CData(2:end, :) = repmat(compColor, length(compTimeValues)-1, 1); % RALMPC (Computation Time) in orange

ylabel('Computation Time [s]', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k'); % Y-label in black
set(gca, 'YColor', 'k'); % Left Y-axis color set to black

% ---- Plot Cost on Right Y-Axis ----
yyaxis right;
b2 = bar(xPos + 0.2, costValues, 0.4, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2.5); % Thicker line for cost
b2.CData(1, :) = refColor; % RAMPC (Cost) in blue
b2.CData(2:end, :) = repmat(compColor, length(costValues)-1, 1); % RALMPC (Cost) in orange

ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k'); % Y-label in black
set(gca, 'YColor', 'k'); % Right Y-axis color set to black
ylim([100, 160]); % Limit the right Y-axis

% ---- Additional Formatting ----
% Set x-axis labels
xticks(xPos);
xticklabels(categories);


% Add a horizontal line for the optimal cost (colorblind-friendly)
yline(J_OS, '--', 'Optimal Cost', 'Color', '#0072B2', ...
    'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);

% Grid and Axis formatting
grid on;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% ---- Create Corrected Legend (Bottom Left) ----
hold on;
h1 = bar(nan, nan, 'FaceColor', refColor, 'EdgeColor', 'none'); % RAMPC (Time)
h2 = bar(nan, nan, 'FaceColor', compColor, 'EdgeColor', 'none'); % RALMPC (Time)
h3 = bar(nan, nan, 'FaceColor', refColor, 'EdgeColor', 'k', 'LineWidth', 2.5); % RAMPC (Cost)
h4 = bar(nan, nan, 'FaceColor', compColor, 'EdgeColor', 'k', 'LineWidth', 2.5); % RALMPC (Cost)
h5 = plot(nan, nan, '--', 'Color', '#0072B2', 'LineWidth', 2); % Optimal Cost Line

legend([h1, h2, h3, h4, h5], ...
    {'RAMPC (Time)', 'RALMPC (Time)', 'RAMPC (Cost)', 'RALMPC (Cost)', 'Optimal Cost'}, ...
    'Location', 'southwest', 'Interpreter', 'latex');

% ---- Export as PDF and EPS ----
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Combined.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Combined.eps');

disp('PaperPlot_Combined.pdf and PaperPlot_Combined.eps');
%% Figure
% Improved plot for CDC - Theta set (Colorblind Friendly)

% Define colorblind-friendly colors from the CUD palette
color_HC0 = [86, 180, 233] / 255;   % Blue
color_RALMPC = [230, 159, 0] / 255; % Orange
color_RAMPC = [0, 158, 115] / 255;  % Green

% Plot the three datasets
figure;
hold on;

% Plot Theta_0^{HC} (first dataset) - Blue
h1 = plot(Theta_HC0, 'LineWidth', 2, 'Color', color_HC0, 'LineStyle', '-', 'DisplayName', 'Theta_0^{HC}');

% Plot Theta_{RALMPC} (last element of RALMPC) - Orange (Solid line)
h2 = plot(Theta_HC_RALMPC{end}{end}, 'LineWidth', 2, 'Color', color_RALMPC, 'LineStyle', '-', 'DisplayName', 'Theta_{RALMPC}');

% Plot Theta_{RAMPC} (last element of RAMPC) - Green (Solid line)
h3 = plot(Theta_HC_RAMPC{end}{end}, 'LineWidth', 2, 'Color', color_RAMPC, 'LineStyle', '-', 'DisplayName', 'Theta_{RAMPC}');

% Add title and labels
xlabel('$\theta_1$', 'Interpreter', 'latex', 'FontSize', 14);  % X-axis label: $\theta_1$
ylabel('$\theta_2$', 'Interpreter', 'latex', 'FontSize', 14);  % Y-axis label: $\theta_2$

% Add grid
grid on;

% Improve axis settings
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Create a manual legend using the line handles h1, h2, h3
legend([h1, h2, h3], ...
    {'$\Theta_0^{HC}$', '$\Theta_{RALMPC}^{HC}$', '$\Theta_{RAMPC}^{HC}$'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);

% Hold off to finish the plot
hold off;

% ---- Export as PDF and EPS ----
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Theta.pdf');
print('-depsc', '/Users/hannes/Documents/MATLAB/Master_Thesis/Linear_RALMPC/RALMPC_LTI_Query/Scenario 2/Plots/PaperPlot_Theta.eps');

disp('PaperPlot_Combined.pdf and PaperPlot_Combined.eps');