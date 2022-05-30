%---------------------------------------------------------------------------
% Free ends sail numerical solution with varying A parameter % 
%---------------------------------------------------------------------------

close all

%---------------------------------------------------------------------------
% Setting up the differential equations %
%---------------------------------------------------------------------------

% Declaring the variables and parameters
syms A C l L s x(s) y(s) Y 

% Declaring the derivatives
Dx = diff(x); % dx/ds
D2x = diff(x,2); % d^2x/ds^2
Dy = diff(y); % dy/ds
D2y = diff(y,2); % d^2y/ds^2

% Position equations of the sail
% d^2x/ds^2
Eq1 = D2x == ...
    A * exp(-x(s)^2/(2*l^2)) * Dy*Dx^2/(Dx^2 + Dy^2)^(3/2) + C*x(s);

% d^2y/ds^2 
Eq2 = D2y == ...
    -A * exp(-x(s)^2/(2*l^2)) * Dx^3/(Dx^2 + Dy^2)^(3/2);

% Converting the second order equations to a system of linear ODEs
[VF,Subs] = odeToVectorField(Eq1, Eq2);
ftotal = matlabFunction(VF,'Vars',{A,C,l,L,s,Y});
% Output: [Y(1), Y(2), Y(3), Y(4)] =  [y, dy/ds, x, dx/ds]




%---------------------------------------------------------------------------
% Parameters %
%---------------------------------------------------------------------------

% All lengths have arbitrary units

L0 = 10; % rope rest length
l = 1e-1*L0; % approximate intensity beam half-width
L = L0/2; % s-parameter maximum, s in [-L, L]
steps = 10001; % Number of integration steps

% A values
A_min = 1e-3;
A_max = 1e1;
logA_min = log10(A_min);
logA_max = log10(A_max);
A_vals = logspace(logA_min,logA_max,7);

% C value - fixed
C = 1e-3;

% Initial conditions at the centre
y_init = 5; % arbitrary maximum y value y(0) = y_init
Dy_init = 0; % zero vertical component, y'(0) = Dy_init
x_init = 0; % arbitrary centre rope x value, x(0) = x_init
Dx_init = 1; % horizontal derivative, x'(0) = Dx_init
ic = [y_init;Dy_init;x_init;Dx_init]; % Initial condition [y,y',x,x']

% Parameter range - right half of the sail
s_range_right = linspace(0,L,steps); 

% Parameter range - left half of the sail
s_range_left = -linspace(0,L,steps); 

% Integrate left and right, outwards from the centre of the sail such that 
% we can use the centre of the rope as the initial condition.

% Solver options for error tolerance
rel_tol = 1e-8;
abs_tol = 1e-10;
options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);




%---------------------------------------------------------------------------
% Solving %
%---------------------------------------------------------------------------

% Plotting parameters
line_width = 2;
axes_fontsize = 24;
colourbar_title_fontsize = 18;
x_range = 6;

set(0,'DefaultAxesFontSize',axes_fontsize)
set(0,'defaultAxesFontName','Times')


% Setting up the figure for plotting
f = figure(1);
    f.Color = "white"; hold on
colour_shift = 40;
c = jet(length(A_vals)+colour_shift);
orange_colour_shift = 29;


% Solving for different A values
for j = 1:length(A_vals)
    A = A_vals(j);
    
    % Integrating over the right half
    [~,Y_right] = ode78(@(s,Y) ftotal(A,C,l,L,s,Y), ...
        s_range_right, ic, options);
    x_right = Y_right(:,3).';
    y_right = Y_right(:,1).';
    % Integrating over the left half
    [~,Y_left] = ode78(@(s,Y) ftotal(A,C,l,L,s,Y), ...
        s_range_left, ic, options);
    x_left = Y_left(:,3).';
    y_left = Y_left(:,1).';
    
    % Full range in the correct order [-L,L]. Remove the doubled entry at 
    % x = x_init, y = y_init in the middle.
    x = [flip(x_left) x_right(2:end)];
    y = [flip(y_left) y_right(2:end)];
    
    % Plotting the equilibrium shape for the current A value
    p = plot(x,y,'Color',c(j+orange_colour_shift,:), ...
        'LineWidth',line_width);
    
    % Limiting the plot axes to the widest sail
    if j == 1
        xlim([x(1),x(end)])
    end
end


% Axes labelling
formatSpec = 'C = %.0e; L_{0} = %.0f; L = %g; l/L_{0} = %g' ;
A_sub = sprintf(formatSpec,C,L0,L,l/L0);
title(A_sub,'FontWeight','Normal')    
xlabel("$x$ (arb. units)",'interpreter','latex')
ylabel("$y$ (arb. units)",'interpreter','latex')


% Gaussian beam profile
beam_range = linspace(-x_range,x_range,200);
intensity_profile = exp(-beam_range.^2./(2*l^2));
intensity_plot = plot(beam_range, max(y)*intensity_profile, ':k', ...
    'Linewidth', line_width);
xlim([-x_range,x_range])    
leg = legend(intensity_plot,'Beam profile');
    leg.FontSize = 20;
hold off

formatSpec = '%.0e' ;
A_min_leg = sprintf(formatSpec,A_min);
A_max_leg = sprintf(formatSpec,A_max);

c_orange = c(orange_colour_shift+1:length(A_vals)+orange_colour_shift, :);
colormap(c_orange)
colour_grad = colorbar('Ticks',[0,1], ...
         'TickLabels',{A_min_leg,A_max_leg}, ...
         'Direction','Reverse');
title(colour_grad, 'A', 'FontSize', colourbar_title_fontsize)