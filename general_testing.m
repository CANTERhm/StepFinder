%% curve data

clearvars -except handles main

table = handles.guiprops.Features.edit_curve_table;
curvename = table.UserData.CurrentCurveName;
curve_data = handles.curveprops.(curvename).RawData.CurveData;
clamp_x = curve_data.Segment5.time;
clamp_y = curve_data.Segment5.vDeflection;

%% steps
finder = StepFinder(clamp_x, clamp_y);
finder.window_width = 500;
finder.smoothing_sigma = 3;
finder.peak_threshold = 0.5;
finder.step_refinement = 1;
finder = finder.SmoothData();
finder = finder.StepSearch();

%% plot theta and data
theta = finder.step_results.theta;
indices = finder.step_results.indices;
w = round(finder.window_width/2);
t_x = [];
t_y = [];
t_x_2 = [];
t_y_2 = [];
for i = 1:length(finder.step_indices)
    t_x = [t_x clamp_x(finder.step_indices(i))];
    t_y = [t_y clamp_y(finder.step_indices(i))];
end

finder = finder.RecalculateStep();

for i = 1:length(finder.step_indices)
    t_x_2 = [t_x_2 clamp_x(finder.step_indices(i))];
    t_y_2 = [t_y_2 clamp_y(finder.step_indices(i))];
end

figure()
hold on
plot(clamp_x, clamp_y);
plot(clamp_x, finder.y_conv);
scatter(t_x, t_y, 'Marker', 'o',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'red',...
    'SizeData', 20);

scatter(t_x_2, t_y_2, 'Marker', 'o',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'green',...
    'SizeData', 20);

hold off
grid on
grid minor
plottools