%---------------------------------------------------------------------------
% Free ends sail - checking error in the numerical solution % 
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

% Initial conditions at the centre
y_init = 10; % arbitrary maximum y value y(0) = y_init
Dy_init = 0; % zero vertical component, y'(0) = Dy_init
x_init = 0; % arbitrary centre rope x value, x(0) = x_init
Dx_init = 1; % horizontal derivative, x'(0) = Dx_init
ic = [y_init;Dy_init;x_init;Dx_init]; % Initial condition [y,y',x,x']
steps = 10001; % Number of integration steps

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

% A and C values
A = 1;
C = 1e-3;

% Right half
[~,Y_right] = ode78(@(s,Y) ftotal(A,C,l,L,s,Y), ...
    s_range_right, ic, options);
x_right = Y_right(:,3).';
Dx_right = Y_right(:,4).';
y_right = Y_right(:,1).';
Dy_right = Y_right(:,2).';
% Left half
[~,Y_left] = ode78(@(s,Y) ftotal(A,C,l,L,s,Y), ...
    s_range_left, ic, options);
x_left = Y_left(:,3).';
Dx_left = Y_left(:,4).';
y_left = Y_left(:,1).';
Dy_left = Y_left(:,2).';

% Full range in the correct order [-L,L]. Remove the doubled entry at 
% x = x_init, y = y_init in the middle.
s = [flip(s_range_left) s_range_right(2:end)]; % domain for x(s) and y(s)
x = [flip(x_left) x_right(2:end)];
y = [flip(y_left) y_right(2:end)];
Dx = [flip(Dx_left) Dx_right(2:end)];
Dy = [flip(Dy_left) Dy_right(2:end)];
Dx_vec_length = sqrt(Dx.^2 + Dy.^2); % length of the tangent vector




%------------------------------------------------------------------------------------------
% Numerical derivatives/solution checking %
%------------------------------------------------------------------------------------------

% Checking that the solution's numerical derivatives solve the original
% differential equations. Take the difference between x_ss, y_ss and the 
% RHS of the sail equations. Expect to see absolute error with lower 
% magnitude than the terms of the equations being solved. 

% Numerical derivatives
s_1 = s(1:end-1); % s truncated by 1 to match the first difference
Dx2 = diff(Dx)./diff(s); % d^2x/ds^2 numerical
Dy2 = diff(Dy)./diff(s); % d^2y/ds^2 numerical
mean_Dx2 = mean(abs(Dx2)); % mean value across the second x-derivative
mean_Dy2 = mean(abs(Dy2)); % mean value across the second y-derivative


% Substituting numerical values into the RHS of the sail equations
x_solver_RHS = A*exp(-x(1:end-1).^2/(2*l^2)) ...
    .* Dy(1:end-1).*Dx(1:end-1).^2 ...
    ./ (Dx(1:end-1).^2+Dy(1:end-1).^2).^(3/2) ...
    + C*x(1:end-1); 
y_solver_RHS = -A*exp(-x(1:end-1).^2/(2*l^2)) ...
    .* Dx(1:end-1).^3 ...
    ./ (Dx(1:end-1).^2+Dy(1:end-1).^2).^(3/2);

% Absolute errors at each s value
x_error = Dx2 - x_solver_RHS;
y_error = Dy2 - y_solver_RHS;


% Mean errors across all s values
mean_x_error = mean(abs(x_error));
mean_y_error = mean(abs(y_error));
mean_x_perc_error = 100*mean_x_error/mean_Dx2;
mean_y_perc_error = 100*mean_y_error/mean_Dy2;
disp(['Mean percentage error in x: ',num2str(mean_x_perc_error),'%.'])
disp(['Mean percentage error in y: ',num2str(mean_y_perc_error),'%.'])


% Plotting the second derivative and absolute error
% Setting up the figure for plotting
f = figure(3);
error_graphs = tiledlayout(2,1);
error_graphs.Padding = 'compact';
error_graphs.TileSpacing = 'compact';
xlabel(error_graphs,"$s$",'interpreter','latex','FontSize',axes_fontsize)

% Plotting the second derivative x_ss and y_ss against s
nexttile
    f.Color = "white";
plot(s_1,Dx2,'-b','LineWidth',2)    
    ylabel("Second derivative")
    set(gca,'box','off')
    hold on
plot(s_1,Dy2,'-r','LineWidth',2)
    legend("$x''(s)$","$y''(s)$",'interpreter','latex')
    hold off
dim = [.2 .5 .3 .3];
str = '(a)';
text(4,-0.75,str, 'FontName','Times', 'FontSize',24);

% Plotting the absolute error against s
nexttile
    f.Color = "white"; hold on
plot(s_1,x_error,'-r','LineWidth',2)
    ylabel("Absolute error (arb. units)")
    hold on
plot(s_1,y_error,'-b','LineWidth',2)
    hold off
str = '(b)';
text(4,-1.5e-4,str, 'FontName','Times', 'FontSize',24);





%------------------------------------------------------------------------------------------
% Length of the rope %
%------------------------------------------------------------------------------------------

% Trapezoidal rule integration of the function sqrt{x'^2 + y'^2} over the
% s-range implicit in the values Dx and Dy. From a definition in 
% mathematics, gives the arc-length.
rope_lengths = cumtrapz(s,Dx_vec_length);
final_rope_length = rope_lengths(end);
disp(['Final rope length: ',num2str(final_rope_length),' arbitrary units.'])