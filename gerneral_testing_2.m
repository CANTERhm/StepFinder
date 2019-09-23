%% curve data

clearvars -except handles main

table = handles.guiprops.Features.edit_curve_table;
curvename = table.UserData.CurrentCurveName;
curve_data = handles.curveprops.(curvename).RawData.CurveData;
clamp_x = curve_data.Segment5.time;
clamp_y = curve_data.Segment5.vDeflection;

%% steps
finder = StepFinder(clamp_x, clamp_y);
finder.window_width = 100;
finder.smoothing_sigma = 3;
finder.peak_threshold = 0.5;
finder.step_refinement = 1;
finder = finder.SmoothData();
finder = finder.StepSearch();
finder = finder.RecalculateStep();

%%
indices = finder.recalculate_step.indices;
g = finder.recalculate_step.g;
f = finder.recalculate_step.f;
pos = finder.recalculate_step.pos;
theta = finder.recalculate_step.theta;


figure()
hold on
plot(clamp_x, clamp_y);
plot(clamp_x, finder.y_conv);
grid on
grid minor
plottools

ind = indices(:,1);
xdat = clamp_x(ind(1):ind(2));

lb_x = [clamp_x(ind(1)) clamp_x(ind(1))];
lb_y = [min(clamp_y) max(clamp_y)];

rb_x = [clamp_x(ind(2)) clamp_x(ind(2))];
rb_y = [min(clamp_y) max(clamp_y)];

mid_x = [clamp_x(pos(1)) clamp_x(pos(1))];
mid_y = [min(clamp_y) max(clamp_y)];

f_i = f(:,1);

g_i = g(:,1);

lb_obj = plot(lb_x, lb_y, 'k-');
rb_obj = plot(rb_x, rb_y, 'k-');
mid_obj = plot(mid_x, mid_y, 'k--');
f_obj = plot(f_i, xdat, 'r-');
g_obj = plot(g_i, xdat, 'b-');

for i = 1:length(indices)
    ind = indices(:,i);
    f_i = f(:,i);
    g_i = g(:,i);
    mid_x = [clamp_x(pos(i)) clamp_x(pos(i))];
    lb_x = [clamp_x(ind(1)) clamp_x(ind(1))];
    rb_x = [clamp_x(ind(2)) clamp_x(ind(2))];
    xdat = clamp_x(ind(1):ind(2));
    
    lb_obj.XData = lb_x;
    rb_obj.XData = rb_x;
    mid_obj.XData = mid_x;
    f_obj.XData = xdat;
    f_obj.YData = f_i;
    g_obj.XData = xdat;
    g_obj.YData = g_i;
    
    fprintf('\ntheta: %g \t index: %d\n', theta(i), i);
    
end
