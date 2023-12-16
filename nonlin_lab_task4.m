%% Nonlinear Control Lab Task 4 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Define the system equations
f1 = x2;
f2 = -(B_l/J_l)*x2 - (k/J_l)*(x1 - x3) - (m*g*l/J_l)*cos(x1);
f3 = x4;
f4 = (k/J_m)*(x1 - x3) - (B_m/J_m)*x4 + nu/J_m;

% Define the nonlinear system
f = [f1; f2; f3; f4];

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
det_A_lin = det(A_lin)

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
% Ensure you have the Symbolic Math Toolbox installed and active
% Define the symbolic variables
syms g l m J_l J_m B_l B_m k real

% Define the system matrices A and B
A = [0, 1, 0, 0;
     sqrt(2)*g*l*m/(2*J_l) - k/J_l, -B_l/J_l, k/J_l, 0;
     0, 0, 0, 1;
     k/J_m, 0, -k/J_m, -B_m/J_m];

B = [0; 0; 0; 1/J_m];

% Define the number of iterations, which should be n for an n-dimensional system
% to include the nth iterated Lie bracket
n = 4; % Assuming the system is 4-dimensional

% Preallocate cell array to store the vector fields
vec_fields = cell(1, n);

% Calculate the vector fields ad_a^0 b, ..., ad_a^(n-1) b
for k = 0:n-1  % Updated range to include k = 3
    % Compute ad_a^k b using the simplified formula for linear systems
    ad_a_b = (-1)^k * A^k * B;
    
    % Store the result in the cell array
    vec_fields{k+1} = simplify(ad_a_b);
end

% Display the results
for i = 1:length(vec_fields)
    fprintf('ad_a^%d b:\n', i-1);
    disp(vec_fields{i});
end

% Create a matrix with the vector fields as columns
V = [vec_fields{:}];  % Concatenate the cell array contents into a matrix

% Calculate the rank of the matrix V
rank_V = rank(V);

% Check if the vectors are linearly independent
if rank_V == size(V, 2)
    fprintf('Constraint 1 = satisfied\n');
else
    fprintf('The vectors are not linearly independent.\n');
end

%% Check whether the controllability and involutivity conditions are satisfied in 𝐷 containing x^0

is_involutory = true;
for i = 1:length(vec_fields)-1  % Only need to go up to n-2
    for j = i+1:length(vec_fields)-1  % Only need to go up to n-2
        % Compute the Lie bracket of the i-th and j-th vector fields
        lie_bracket = jacobian(vec_fields{j}, x) * vec_fields{i} - jacobian(vec_fields{i}, x) * vec_fields{j};
        
        % Check if this Lie bracket is in the span of the original set
        % Simplified check, assumes the rank is full and the set is already a basis
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

% Define the state vector
x = [theta_l; theta_l_dot; theta_m; theta_m_dot];

% Define the equilibrium point x^0
x0 = [pi/4; 0; 0; 0];

% Define the vector fields a(x) and b(x)
A = [0, 1, 0, 0;
     sqrt(2)*g*l*m/(2*J_l) - k/J_l, -B_l/J_l, k/J_l, 0;
     0, 0, 0, 1;
     k/J_m, 0, -k/J_m, -B_m/J_m];
B = [0; 0; 0; 1/J_m];

% Compute the Lie derivatives L_b L_a^i c for i = 0 to n-2
n = 4; % Assuming the system is 4-dimensional
L_b_L_a_c = c; % Initialize the Lie derivative of c with respect to a

for i = 0:n-2
    L_b_L_a_c = jacobian(L_b_L_a_c, x) * B; % Compute L_b (L_a^i c)
    % Evaluate the Lie derivative at the equilibrium point x^0
    L_b_L_a_c_at_x0 = double(subs(L_b_L_a_c, x, x0));
    if L_b_L_a_c_at_x0 == 0
        fprintf('L_b L_a^%d c = 0 at x^0\n', i);
    else
        fprintf('L_b L_a^%d c is not zero at x^0\n', i);
        break; % If any of the first n-2 Lie derivatives are non-zero, the system is not IO linearizable
    end
    
    if i < n-2
        L_a_c = jacobian(L_b_L_a_c, x) * A; % Compute the next L_a^i c for the next iteration
    end
end

% Compute the (n-1)th Lie derivative L_b L_a^(n-1) c
L_b_L_a_n_minus_1_c = jacobian(L_a_c, x) * A * B;
L_b_L_a_n_minus_1_c_at_x0 = double(subs(L_b_L_a_n_minus_1_c, x, x0));

% Check if L_b L_a^(n-1) c is not zero at x^0
if L_b_L_a_n_minus_1_c_at_x0 ~= 0
    fprintf('L_b L_a^(n-1) c is not zero at x^0, so the system is input-output linearizable at this point.\n');
else
    fprintf('L_b L_a^(n-1) c = 0 at x^0, so the system is not input-output linearizable at this point.\n');
end

