%% Nonlinear Control Lab
syms theta_l theta_l_dot theta_m theta_m_dot nu m l g J_l J_m B_l B_m k real

% Define the state variables
x1 = theta_l;       % theta_l
x2 = theta_l_dot;   % dot(theta_l)
x3 = theta_m;       % theta_m
x4 = theta_m_dot;   % dot(theta_m)

% Define the system equations
f1 = x2;
f2 = -(B_l/J_l)*x2 - (k/J_l)*(x1 - x3) - (m*g*l/J_l)*cos(x1);
f3 = x4;
f4 = (k/J_m)*(x1 - x3) - (B_m/J_m)*x4 + nu/J_m;

% Define the nonlinear system
f = [f1; f2; f3; f4];

% Define the point of linearization
theta_l_eq = pi/4;
x_eq = [theta_l_eq; 0; theta_l_eq; 0];  % Equilibrium point
u_eq = 0;  % Assuming no input at the equilibrium

% Substitute parameters (use your own values)
param_values = [m, l, g, J_l, J_m, B_l, B_m, k];
param_numeric = [0.3, 0.3, 9.8, 4e-4, 4e-4, 0.0, 0.015, 0.8]; % Example values

% Correctly combine all variables and values for substitution
all_variables = [x1, x2, x3, x4, nu, param_values];
all_values = [x_eq', u_eq, param_numeric]; % Ensure all_values is a row vector

% Linearize the system
A = double(subs(jacobian(f, [x1, x2, x3, x4]), all_variables, all_values));
B = double(subs(jacobian(f, nu), all_variables, all_values));

%K = eye(4);

% Example desired pole locations (adjust these based on your system)
poles = [-1 -1.1 -1.2 -1.3];

% Compute the feedback gain matrix
K = place(A, B, poles);

x_init = [0.1; 0.2; 0.3; 0.4];  % New initial values
x = x_init;  % Use these initial values for x

x_ref = [pi/4; 0 ; 0; 0];

% Display the linearized system matrices
disp('Linearized system matrix A:');
disp(A);
disp('Linearized system matrix B:');
disp(B);