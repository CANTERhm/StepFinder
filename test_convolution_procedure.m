%% curve data

clearvars -except handles main

table = handles.guiprops.Features.edit_curve_table;
curvename = table.UserData.CurrentCurveName;
curve_data = handles.curveprops.(curvename).RawData.CurveData;
clamp_x = curve_data.Segment5.time;
clamp_y = curve_data.Segment5.vDeflection;


%% plot for smoothed data
fig = figure();
ax = axes(fig, 'NextPlot', 'add');
plot(ax, clamp_x, clamp_y, 'DisplayName', 'Clamp-Data');
grid on
grid minor
plottools;

%% smooth data
finder = StepFinder(clamp_x, clamp_y);
finder = finder.SmoothData();
plot(ax, finder.x_data, finder.y_conv, 'r-', 'DisplayName', 'Smoothded Data');