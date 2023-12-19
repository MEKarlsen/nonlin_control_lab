%% Nonlinear Control Lab Task 4 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc 
%clear

%% Init

syms theta_l theta_l_dot theta_m theta_m_dot nu m l g J_l J_m B_l B_m k real

param_values = [m, l, g, J_l, J_m, B_l, B_m, k];
param_numeric = [0.3, 0.3, 9.8, 4e-4, 4e-4, 0.0, 0.015, 0.8]; 

% Define the state variables
x1 = theta_l;       % theta_l
x2 = theta_l_dot;   % dot(theta_l)
x3 = theta_m;       % theta_m
x4 = theta_m_dot;   % dot(theta_m)

x = [x1; x2; x3; x4];

% Define the system equations
f1 = x2;
f2 = -(B_l/J_l)*x2 - (k/J_l)*(x1 - x3) - (m*g*l/J_l)*cos(x1);
f3 = x4;
f4 = (k/J_m)*(x1 - x3) - (B_m/J_m)*x4 + nu/J_m;

% Define the nonlinear system
f = [f1; f2; f3; f4];

% Define the system equations
a1 = x2;    
a2 = -(B_l/J_l)*x2 - (k/J_l)*(x1 - x3) - (m*g*l/J_l)*cos(x1);
a3 = x4;
a4 = (k/J_m)*(x1 - x3) - (B_m/J_m)*x4;

% Define the nonlinear system
a = [a1; a2; a3; a4];
b = [0; 0; 0; 1/J_m];

% Define the point of linearization
theta_l_eq = pi/4;
theta_m_eq = pi/2;
x_eq = [theta_l_eq; 0; 0; 0]; 
%u_eq = m*g*l*cos(theta_l_eq);
u_eq = - k * (theta_l_eq - theta_m_eq);

% Correctly combine all variables and values for substitution
all_variables = [x1, x2, x3, x4, nu, param_values];
all_values = [x_eq', u_eq, param_numeric]; % Ensure all_values is a row vector

% Linearize the system
A_lin = double(subs(jacobian(f, [x1, x2, x3, x4]), all_variables, all_values));
B_lin = double(subs(jacobian(f, nu), all_variables, all_values));

% Calculate the determinant of A_lin
det_A_lin = det(A_lin);

% Check if the matrix is singular
if det_A_lin == 0
    disp('A_lin is singular');
else
    disp('A_lin is non-singular');
end

x_init = [0 0 theta_m_eq 0]';  
x_ref = [theta_l_eq 0 theta_m_eq 0]';

%%% Regulators %%%
%%%%%%%%%%%%%%%%%%

%% Full state feedback controller

% Example desired pole locations (adjust these based on your system)
poles = [-2000 -2500 -3000 -3500];

% Compute the feedback gain matrix
K_poles = place(A_lin, B_lin, poles);

%% LQR w/ augmented integral state for Theta_l 
% Original state-space model
A = A_lin;  % Your original A matrix
B = B_lin;  % Your original B matrix
C = [1 0 0 0];  % Assumes everything is measured
D = 0;  % Assume no direct feedthrough for simplicity

% Augment the state-space model to include the integral of the error
A_aug = [A, zeros(4, 1); -C, 0];
B_aug = [B; 0];

% Choose new Q and R matrices for the augmented system
Q_aug = diag([10, 200, 1, 50, 500]);  % The last element is for the integrator state
R_aug = 100;  % Control-effort weighting

% Compute the LQR gain for the augmented system
[K_aug, S_aug, e_aug] = lqr(A_aug, B_aug, Q_aug, R_aug);

% Extract the feedback gain and the integrator gain from the augmented gain matrix
K_2 = K_aug(1:4);  % Feedback gain for the original states
K_i = 1*abs(K_aug(5));  % Integrator gain

%% Define the output
y = pi/4 - theta_l;

%% Compute the Lie Derivatives
L_f_y = diff(y, theta_l) * f1 + diff(y, theta_l_dot) * f2;
L_g_L_f_y = diff(L_f_y, theta_l) * f1 + diff(L_f_y, theta_l_dot) * f2 + diff(L_f_y, theta_m) * f3 + diff(L_f_y, theta_m_dot) * f4 + diff(L_f_y, nu);

%% Determine the Relative Degree
% The relative degree is the number of Lie derivatives needed until the input appears
% In this case, check if L_g_L_f_y is non-zero (which means input appears)

if L_g_L_f_y == 0
    disp('System is not I/O linearizable at this point.');
else
    disp('System is I/O linearizable at this point.');
end

%% Evaluate the expressions at the equilibrium point
L_f_y_at_eq = double(subs(L_f_y, all_variables, all_values));
L_g_L_f_y_at_eq = double(subs(L_g_L_f_y, all_variables, all_values));

disp('L_f_y at equilibrium:');
disp(L_f_y_at_eq);
disp('L_g_L_f_y at equilibrium:');
disp(L_g_L_f_y_at_eq);

%% Construct the vector fields {b, ad_a^1 b  ...} for the given regular system 
% Assuming the system is 4-dimensional
n = 4;

% Preallocate cell array to store the vector fields
vec_fields = cell(1, n);

% Compute the vector fields ad_a^0 b, ..., ad_a^(n-1) b
vec_fields{1} = b; % ad_a^0 b is just b
for k = 1:n-1
    % Compute ad_a^k b using the recursive definition
    vec_fields{k+1} = simplify(jacobian(vec_fields{k}, x) * a - jacobian(a, x) * vec_fields{k});
end

% Display the results
for i = 1:length(vec_fields)
    fprintf('ad_a^%d b:\n', i-1);
    disp(vec_fields{i});
end

% Create a matrix with the vector fields as columns
V = [vec_fields{:}];

% Calculate the rank of the matrix V
rank_V = rank(V);

% Check if the vectors are linearly independent
if rank_V == size(V, 2)
    fprintf('Constraint 1 = satisfied\n');
else
    fprintf('The vectors are not linearly independent.\n');
end

%% Check whether the controllability and involutivity conditions are satisfied in 𝐷 containing x^0

% Check involutivity
is_involutory = true;
for i = 1:n-2
    for j = i+1:n-2
        % Compute the Lie bracket of the i-th and j-th vector fields
        lie_bracket = jacobian(vec_fields{j}, x) * vec_fields{i} - jacobian(vec_fields{i}, x) * vec_fields{j};
        
        % Check if this Lie bracket is in the span of the original set
        coeff = V\lie_bracket;
        if norm(V*coeff - lie_bracket) > 1e-6 % A tolerance for numerical error
            is_involutory = false;
            break;
        end
    end
    if ~is_involutory
        break;
    end
end

if is_involutory
    fprintf('The set of vector fields is involutive in D.\n');
else
    fprintf('The set of vector fields is not involutive in D.\n');
end

%% Check output function
% Define the output function c(x)
c = pi/4 - theta_l;

% Define the equilibrium point x^0
x0 = [pi/4; 0; 0; 0];

% First Lie derivatives
Lb_c = lie_derivative(b, c, x);
La_c = lie_derivative(a, c, x);

% Higher-order Lie derivatives
Lb_La_c = lie_derivative(b, La_c, x);
La_La_c = lie_derivative(a, La_c, x);
Lb_La_La_c = lie_derivative(b, La_La_c, x);
La_La_La_c = lie_derivative(a, La_La_c, x);
Lb_La_La_La_c = lie_derivative(b, La_La_La_c, x);
La_La_La_La_c = lie_derivative(a, La_La_La_c, x);

% Relative degree
r = 4; % Lb_La^(r-1)_c = Lb_La^(4-1)_c = -k/(J_l*J_m) =/= 0, 

disp('La^k derivatives:');
disp(La_c);
disp(La_La_c);
disp(La_La_La_c);
disp(La_La_La_La_c);

% Evaluate at the equilibrium point x0 if necessary
Lb_c_at_x0 = subs(Lb_c, x, x0);
La_c_at_x0 = subs(La_c, x, x0);
Lb_La_c_at_x0 = subs(Lb_La_c, x, x0);
Lb_La_La_c_at_x0 = subs(Lb_La_La_c, x, x0);
Lb_La_La_La_c_at_x0 = subs(Lb_La_La_La_c, x, x0);
La_La_La_La_c_at_x0 = subs(La_La_La_La_c, x, x0);

% Substitute the numerical values into the expression
La_La_La_La_c_at_x0_numeric = double(subs(La_La_La_La_c_at_x0, param_values, param_numeric));
Lb_La_La_La_c_at_x0_numeric = double(subs(Lb_La_La_La_c_at_x0, param_values, param_numeric));

% Display the evaluated Lie derivatives at the equilibrium point
%disp('Lie derivatives evaluated at the equilibrium point:');
%disp(Lb_c_at_x0);
%disp(La_c_at_x0);
%disp(Lb_La_c_at_x0);
%disp(Lb_La_La_c_at_x0);
%disp(Lb_La_La_La_c_at_x0);

syms u; % Control input

% Calculate the Lie derivatives for k = 1 to 4
for k = 1:4
    L_a_k_c = higher_order_lie_derivative(a, c, x, k);
    if k > 1
        L_b_L_a_k_minus_1_c = lie_derivative(b, higher_order_lie_derivative(a, c, x, k-1), x);
    else
        L_b_L_a_k_minus_1_c = 0; % L_b_L_a^0_c is zero since L_a^0_c = c and L_b_c has already been calculated if needed
    end
    y_k = L_a_k_c + u * L_b_L_a_k_minus_1_c;
    
    % Evaluate at the equilibrium point x0 if necessary
    %y_k_at_x0 = subs(y_k, x, x0);
    
    % Display the k-th derivative
    fprintf('The %d-th derivative of y is:\n', k);
    disp(y_k);
end

v = La_La_La_La_c + u * Lb_La_La_La_c;
u = solve(v, u);

% The new coordinates
x_tilde_1 = c;
x_tilde_2 = La_c;
x_tilde_3 = La_La_c;
x_tilde_4 = La_La_La_c;

% Define the new state vector in terms of the new coordinates
x_tilde = [x_tilde_1; x_tilde_2; x_tilde_3; x_tilde_4];

% Calculate the Jacobian of the new state vector with respect to the original state vector
J_x_tilde = jacobian(x_tilde, x);

% Display the Jacobian matrix
disp('Jacobian of x_tilde with respect to x:');
disp(J_x_tilde);

% Substitute the numeric values at the equilibrium point if necessary
J_x_tilde_at_eq = double(subs(J_x_tilde, all_variables, all_values));

% Display the Jacobian evaluated at the equilibrium point
disp('Jacobian of x_tilde at the equilibrium point:');
disp(J_x_tilde_at_eq);  

% Substitute the numeric values at the equilibrium point if necessary
x_tilde_1_at_eq = double(subs(x_tilde_1, all_variables, all_values));
x_tilde_2_at_eq = double(subs(x_tilde_2, all_variables, all_values));
x_tilde_3_at_eq = double(subs(x_tilde_3, all_variables, all_values));
x_tilde_4_at_eq = double(subs(x_tilde_4, all_variables, all_values));

% You can also perform the check at the equilibrium point
rank_J_x_tilde_at_eq = rank(J_x_tilde_at_eq);

% Check if the Jacobian at the equilibrium point is full rank
if rank_J_x_tilde_at_eq == length(x)
    disp('The Jacobian matrix at the equilibrium point is full rank.');
else
    disp('The Jacobian matrix at the equilibrium point is NOT full rank.');
end

% Display the new coordinates evaluated at the equilibrium point
disp('The new coordinates at the equilibrium point are:');
disp(x_tilde_1_at_eq);
disp(x_tilde_2_at_eq);
disp(x_tilde_3_at_eq);
disp(x_tilde_4_at_eq);

% Calculate the time derivatives of the new state variables (x_tilde)
x_tilde_dot_1 = x_tilde_2;
x_tilde_dot_2 = x_tilde_3;
x_tilde_dot_3 = x_tilde_4;
syms u_old;
%x_tilde_dot_4 = La_La_La_La_c_at_x0 + u_old*Lb_La_La_La_c_at_x0; % use
%equations linearized around qeuilibrium when creating simulink!
x_tilde_dot_4 = La_La_La_La_c + u_old*Lb_La_La_La_c;

x_tilde_dot = [x_tilde_dot_1; x_tilde_dot_2; x_tilde_dot_3; x_tilde_dot_4];
disp('x_tilde_dot: ');
disp(x_tilde_dot);

%% LQR regulator

% Define the state-space matrices in the new coordinates
A_tilde = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0];
B_tilde = [0; 0; 0; 1];
C_tilde = [1 0 0 0];
D_tilde = 0;

% Choose the weighting matrices Q and R
Q = diag([1, 1, 1, 1]);  % Example: unit weights on all states
R = 1;  % Example: unit weight on input effort

% Solve the LQR problem
[K_lin_feedback, S, e] = lqr(A_tilde, B_tilde, Q, R);

% Display the optimal gain matrix
disp('The LQR gain matrix K is:');
disp(K_lin_feedback);

%% Functions

% Define a function to calculate the Lie derivative
function Lf_g = lie_derivative(f, g, x)
    % Calculates the Lie derivative of g with respect to vector field f
    Lf_g = 0;
    for i = 1:length(x)
        Lf_g = Lf_g + diff(g, x(i))*f(i);
    end
end

function Lf_g = higher_order_lie_derivative(f, g, x, order)
    % Computes the higher-order Lie derivative of g with respect to vector field f
    Lf_g = g;
    for k = 1:order
        Lf_g = lie_derivative(f, Lf_g, x);
    end
end
