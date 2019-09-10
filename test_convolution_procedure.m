curve_data = handles.curveprops.FC_171214M2C1_0245.RawData.CurveData;
clamp_x = curve_data.Segment5.seriesTime;
clamp_y = curve_data.Segment5.vDeflection;

fig = figure();
ax = axes(fig, 'NextPlot', 'add');
plot(ax, clamp_x, clamp_y, 'DisplayName', 'Clamp-Data');
grid on
grid minor
legend 

% calculate convolution
% x_gauss = linspace(-1,1,100);
% % g_vec = gauss(x_gauss,0,0.01);
% g_vec = normpdf(x_gauss,0,50);
% y_conv = conv(clamp_y, g_vec, 'same');
y_conv = smoothdata(clamp_y,'gaussian','SmoothingFactor',0.02);

plot(ax, clamp_x, y_conv, 'r-', 'DisplayName', 'Convolved Data');


function g_vec = gauss(x_vec, mu, sig)
%     g_vec = exp(-(x_vec-mu).^2./(2*sig^2));
    g_vec = (1/sqrt(2*pi*sig^2)).*exp(-(x_vec-mu).^2./(2*sig^2));
end
