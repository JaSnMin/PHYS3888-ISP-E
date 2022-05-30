%-----------------------------------------------------------------------
% Fixed ends sail, plotting with varying A parameter %
%-----------------------------------------------------------------------

close all

%-----------------------------------------------------------------------
% Parameters %
%-----------------------------------------------------------------------

% Intensity parameter, A
A_min = 1e-3;
A_max = 1e-1;
logA_min = log10(A_min);
logA_max = log10(A_max);
A_vals = logspace(logA_min,logA_max,7);

L0 = 10; % Rest length of the rope
x_range = L0/2; % Maximum x-value of the sail




%-----------------------------------------------------------------------
% Setting up the figure for plotting %
%-----------------------------------------------------------------------
line_width = 2;
axes_fontsize = 24;
colourbar_title_fontsize = 18;
set(0,'DefaultAxesFontSize',axes_fontsize)
set(0,'defaultAxesFontName','Times')

f = figure(1);
    f.Color = "white"; hold on
colour_shift = 40;
c = jet(length(A_vals)+colour_shift);
orange_colour_shift = 29;




%-----------------------------------------------------------------------
% Plotting the analytic solution across A values
%-----------------------------------------------------------------------
x = linspace(-x_range,x_range,200); % x-values of the sail

% Plotting the sail shape (y) against x for each A value
for j = 1:length(A_vals)
    % Analytic equilibrium sail shape y(x)
    y = (cosh(A_vals(j)*L0/2) - cosh(A_vals(j)*x)); 
    
    % Plotting on the same graph
    p = plot(x,y,'Color',c(j+orange_colour_shift,:), ...
        'LineWidth',line_width);
end


% Plotting the uniform beam profile
intensity_profile = max(y)/2*ones(1,length(x));
intensity_plot = plot(x,intensity_profile,':k', ...
    'LineWidth',line_width);
xlim([-x_range,x_range])    
leg = legend(intensity_plot,'Beam profile');
    leg.FontSize = 20;


% Axes labelling
formatSpec = 'L_{0} = %.0f' ;
A_sub = sprintf(formatSpec,L0);
title(A_sub,'FontWeight','Normal')    
xlabel("$x$ (arb. units)",'interpreter','latex')
ylabel("$y$ (arb. units)",'interpreter','latex')


% Legend labelling
formatSpec = '%.0e' ;
A_min_leg = sprintf(formatSpec,A_min);
A_max_leg = sprintf(formatSpec,A_max);

c_orange = c(orange_colour_shift+1:length(A_vals)+orange_colour_shift, :);
colormap(c_orange)
colour_grad = colorbar('Ticks',[0,1], ...
         'TickLabels',{A_min_leg,A_max_leg}, ...
         'Direction','Reverse');
title(colour_grad, 'A', 'FontSize', colourbar_title_fontsize)