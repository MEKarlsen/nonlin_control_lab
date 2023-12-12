%% Nonlinear Control Lab %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear

%% Init

syms theta_l theta_l_dot theta_m theta_m_dot nu m l g J_l J_m B_l B_m k real

param_values = [m, l, g, J_l, J_m, B_l, B_m, k];
param_numeric = [0.3, 0.3, 9.8, 4e-4, 4e-4, 0.0, 0.015, 0.8]; % Example values

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

x_init = [0 0 theta_m_eq 0]';  

x_ref = [theta_l_eq 0 theta_m_eq 0]';

%%% Regulators %%%
%%%%%%%%%%%%%%%%%%

%% Full state feedback controller

% Example desired pole locations (adjust these based on your system)
poles = [-2 -2.5 -3 -3.5];

% Compute the feedback gain matrix
K = place(A_lin, B_lin, poles);

%% LQR
Q = diag([1000, 0.1, 0.1, 0.1]);  % Define the state-error weighting
R = 2;  % Define the control-effort weighting

% Compute the LQR gain
[K, S, e] = lqr(A_lin, B_lin, Q, R);
% K is the optimal gain matrix

%% LQR w/ augmented integral state for Theta_l 
% Original state-space model
A = A_lin;  % Your original A matrix
B = B_lin;  % Your original B matrix
C = [1 0 1 0];  % Assumes everything is measured
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
K = K_aug(1:4);  % Feedback gain for the original states
K_i = 1*abs(K_aug(5));  % Integrator gain

%% Feedforward Gain (NOTE: useless for this system)
% Compute the pseudo-inverse for feedforward gain
K_ff = pinv(B_lin);

% Compute the feedforward control action
u_ff = K_ff * x_ref;

%% Display the linearized system matrices
disp('Linearized system matrix A:');
disp(A_lin);
disp('Linearized system matrix B:');
disp(B_lin);

%% Notater

% u should be 1x1?
% 2 motorer, 1 til hvert ledd => u burde være 2x1? 

% se på tuning 

% use nonlin sys in simulation and give input to the nonlin system.

% bruk lin model til å lage K.: u = - Kx. Bruker per nå LQR.

% ANTAKELSE: alt kan måles, nevn i prestasjon
% IFT. input; 2 motorer. Se på verdier som er reasonable. Kan argumentere
% med LQR og vekting! 

%Feedforward funker ikke pga B, men lar den være i koden fremdeles.
