clc
clear
close all
%% Load offline computed parameter
load("Parameter_Offline.mat")
%% Declare MPC Parameter
import casadi.*
% Prediction horizon
N_RAMPC=12;
% Number of MPC iterations
mpciterations = 100;
numberitertions=10;
%Option for casadi solver
opts = struct;
opts.ipopt.print_level = 0;  % Suppress solver output
opts.print_time = false;      % Disable timing information
%% Settings for other optimization
options = optimset('Display','none',...
    'TolFun', 1e-10,...
    'MaxIter', 1000,...
    'TolConSQP', 1e-8);
%% Other Parameters
x_0=[0.1;0.1];
theta_star=[1;1];
%For plotting the ellipses
[eig_U, eig_D] = eig(P);
theta_circle = linspace(0, 2*pi, 100);
circle = [cos(theta_circle); sin(theta_circle)];
%Parameter estimation 
M=10;
theta_bar_t{1}=theta_0;
theta_hat_t{1}=theta_0;
Theta_t{1}=Theta_0;
rho_theta_t{1}=rho_theta0;
eta_t{1}=eta_0;
for j=1:M
    Delta{j}=Theta_0;
end
%% Compute initial solution
fprintf("Compute inital solution\n");
[SS_inital,J_wc_inital,X_bar_inital,S_inital,J_inital]=computeinitialsolution(x_0,x_s,u_s,theta_bar_t{end},rho_theta_t{end},eta_t{end},L_B_rho,d_bar,L,l,c,c_xs,P,mpciterations,Q,R,m,c_alpha);
%% Solve the optimal control problem
fprintf("Solve the optimal problem\n");
N_OS=100;
mpciterations_OS=100;
[X_OL_OS,U_OL_OS,J_OS,X_OS,U_OS]=solve_OptimalSolution(x_0,x_s, u_s, Q, R,theta_star,L,l,r,N_OS,mpciterations_OS,m,n);
%% RAMPC without learning
fprintf("Solve the RAMPC problem\n");
theta_bar_RAMPC{1}{1}=theta_0;
theta_hat_RAMPC{1}{1}=theta_0;
Theta_HC_RAMPC{1}{1}=Theta_0;
eta_RAMPC{1}{1}=eta_0;
rho_theta_RAMPC{1}{1}=rho_theta0;
for i=1:numberitertions
    tStart = cputime;
    [X_RAMPC{i},U_RAMPC{i},J_RAMPC{i},X_hat_OL_RAMPC{i},X_bar_OL_RAMPC{i},U_bar_OL_RAMPC{i},S_OL_RAMPC{i},Theta_HC_RAMPC{i},theta_bar_RAMPC{i},theta_hat_RAMPC{i},eta_RAMPC{i},rho_theta_RAMPC{i},t_cpu_RAMPC{i},Delta]=RAMPC(x_0,x_s, u_s, Q, R,P,delta_loc,c_xs,c_alpha,Theta_HC_RAMPC{end}{end},theta_star,theta_bar_RAMPC{end}{end},theta_hat_RAMPC{end}{end},eta_RAMPC{end}{end},rho_theta_RAMPC{end}{end},L_B_rho,d_bar,L,l,Delta,D,options,mu,B_p,p,r,N_RAMPC,mpciterations,m,n,c);
    t_cpu_RAMPC{i} = cputime - tStart;
end
%% Query
fprintf("Query: Solve the RALMPC problem\n");
for j=1:M
    Delta{j}=Theta_0;
end
for N=12:-2:4
    % RALMPC Algorithm 
    SS{1}=SS_inital ;
    Q_func{1}=J_wc_inital;
    fprintf("Iteration:\t%d\n", N);
    [X_RALMPC,U_RALMPC,J_RALMPC,X_hat_OL_RALMPC,X_bar_OL_RALMPC,U_bar_OL_RALMPC,S_OL_RALMPC,J_wc_RALMPC,Theta_HC_RALMPC,theta_bar_RALMPC,theta_hat_RALMPC,eta_RALMPC,t_cpu_RALMPC]=RALMPC(x_0,x_s, u_s, Q, R,P,Theta_0,theta_star,theta_0,eta_0,rho_theta0,L_B_rho,d_bar,L,l,Delta,D,options,mu,B_p,p,r,N,numberitertions,mpciterations,m,n,c,SS,Q_func);
    % Save for plotting
    name="/Users/hannes/Documents/MATLAB/Master_Thesis/Nonlinear_RALMPC/RALMPC Paper example/Data/data_"+"Iter"+string(numberitertions)+"_N"+string(N)+".mat";
    save(name,"N","SS_inital","J_wc_inital","X_bar_inital","S_inital","J_inital","X_RALMPC","U_RALMPC","J_RALMPC","X_hat_OL_RALMPC","X_bar_OL_RALMPC","U_bar_OL_RALMPC","S_OL_RALMPC","J_wc_RALMPC","Theta_HC_RALMPC","theta_bar_RALMPC","eta_RALMPC","t_cpu_RALMPC","X_RAMPC","U_RAMPC","J_RAMPC","X_hat_OL_RAMPC","X_bar_OL_RAMPC","U_bar_OL_RAMPC","S_OL_RAMPC","Theta_HC_RAMPC","theta_bar_RAMPC","eta_RAMPC","t_cpu_RAMPC","X_OL_OS","U_OL_OS","J_OS","X_OS","U_OS")
end