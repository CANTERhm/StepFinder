%% curve data

clearvars -except handles main

table = handles.guiprops.Features.edit_curve_table;
curvename = table.UserData.CurrentCurveName;
curve_data = handles.curveprops.(curvename).RawData.CurveData;
clamp_x = curve_data.Segment5.time;
clamp_y = curve_data.Segment5.vDeflection;

% [clamp_x, clamp_y] = StepSignal(50, 30, 0.5);

%% plot for steps
fig = figure();
ax1 = subplot(2,1,1, 'NextPlot', 'add');
plot(ax1, clamp_x, clamp_y, 'k-', 'DisplayName', 'Clamp-Data');

grid on
grid minor
plottools;

%% calculate steps
finder = StepFinder(clamp_x, clamp_y);
finder.window_width = 50;
finder.smoothing_sigma = 3;
finder.step_refinement = 1;
finder = finder.SmoothData();
finder = finder.StepSearch();
ax2 = subplot(2,1,2);
plot(ax1, clamp_x, finder.y_conv, 'g-', 'DisplayName', 'Smoothed Data');
plot(ax2, clamp_x, finder.theta, 'k-', 'DisplayName', 'theta');
grid on
grid minor

%% plot moving window
indices = finder.step_results.indices;
RSS_g = finder.step_results.RSS_g;
RSS_f = finder.step_results.RSS_f;
f = finder.step_results.f;
g = finder.step_results.g;
% y = finder.step_results.y;
len = length(indices);
width = finder.window_width;
w = round(width/2);
% m = finder.step_results.m;
% m_0 = finder.step_results.m_0;
% t_l = finder.step_results.t_l;
% t_r = finder.step_results.t_r;
% t_0 = finder.step_results.t_0;
% yvec = finder.step_results.y;

a = indices(1,1);
b = indices(2,1);

% width = b-a;
% w = round(width/2);

xfit = clamp_x(a:b);
yfit = f(:,1);
yfit_global = g(:,1);
mid = (w + a);

x_l = xfit(1:round(length(xfit))/2);
x_r = xfit(round(length(xfit))/2+1:end);
x_g = xfit;
f_l = yfit(1:round(length(yfit))/2);
f_r = yfit(round(length(yfit))/2+1:end);
g_i = yfit_global;

% mid = (w + a);
% f_l = f(1:w,1);
% f_r = f(w:end,1);
% x_l = clamp_x(a:w-1+a)';
% x_r = clamp_x(a+w-1:b-1)';
% g_i = g(:,1);
% x_g = clamp_x(a:b-1)';

lb_obj = plot(ax1, [clamp_x(a) clamp_x(a)], [min(clamp_y) max(clamp_y)], 'k-',...
    'DisplayName', 'Left Border');
rb_obj = plot(ax1, [clamp_x(b-1) clamp_x(b-1)], [min(clamp_y) max(clamp_y)], 'k-',...
    'DisplayName', 'Right Border');
mid_obj = plot(ax1, [clamp_x(mid) clamp_x(mid)], [min(clamp_y) max(clamp_y)], 'k--',...
    'DisplayName', 'Evaluation Position');
f_l_obj = plot(ax1, x_l, f_l, 'r-', 'DisplayName', 'Left Piecewise Linear');
f_r_obj = plot(ax1, x_r, f_r, 'r-', 'DisplayName', 'Right Piecewise Linear');
g_obj = plot(ax1, x_g, g_i, 'b-', 'DisplayName', 'Global Windwo Linear');

for i = 430:len
    a = indices(1,i);
    b = indices(2,i);
    xfit = clamp_x(a:b);
    yfit = f(:,i);
    yfit_global = g(:,i);
    mid = (w + a);
    
    x_l = xfit(1:round(length(xfit))/2);
    x_r = xfit(round(length(xfit))/2+1:end);
    x_g = xfit;
    f_l = yfit(1:round(length(yfit))/2);
    f_r = yfit(round(length(yfit))/2+1:end);
    g_i = yfit_global;
%     f_l = f(1:w,i);
%     f_r = f(w:end,i);
%     x_l = clamp_x(a:w-1+a)';
%     x_r = clamp_x(a+w-1:b-1)';
%     g_i = g(:,i);
%     x_g = clamp_x(a:b-1)';
    
    try
        lb_obj.XData = [clamp_x(a) clamp_x(a)];
    catch
    end
    
    try
        rb_obj.XData = [clamp_x(b) clamp_x(b)];
    catch
    end
    
    try
        mid_obj.XData = [clamp_x(mid) clamp_x(mid)];
    catch
    end
    
    try
        f_l_obj.XData = x_l;
        f_l_obj.YData = f_l;
    catch
    end
    
    try
        f_r_obj.XData = x_r;
        f_r_obj.YData = f_r;
    catch
    end
    
    try
        g_obj.XData = x_g;
        g_obj.YData = g_i;
    catch
    end
    
end

function varargout = StepSignal(len, step_pos, step_height)
    x = linspace(0, len, len);
    test_signal = heaviside(x-step_pos)*step_height;
    varargout{1} = x;
    varargout{2} = test_signal;
end