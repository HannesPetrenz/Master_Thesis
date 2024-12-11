clear 
clc
close all
%% Settings
convex=true;
disturbance_deter=true;
%% System Definition
T_s=0.1; %sampling time
%Uncertain Parameter
m=1; %Mass
k=1; %Spring  -->p1
c=0.2; %Damper -->p2
%True Parameter
k_star=0.5;
c_star=0.3;
%correct system matrix
A_star=[1 T_s; -T_s*k_star/m 1-T_s*c_star/m];
B_star = [0; T_s/m];
%System dimensions
n=size(A_star,1);
m=size(B_star,2);
p=2;
%Affine uncertainty description
parametric_unc=0.5; %scaling of the uncertainty s.t. theta_star=[1,-1]
p1_unc = parametric_unc*(k/m);
p2_unc = parametric_unc*(c/m);
A_0 = [1 T_s; -T_s*k/m 1-T_s*c/m];
A_1 = p1_unc*[0 0; T_s 0];
A_2 = p2_unc*[0 0; 0 T_s];
B_0 = [0; T_s/m];
B_1 = zeros(n,m);
B_2 = zeros(n,m);

%Disturbance 
d_max=0.2;
w_max=d_max/m*T_s;%in state space we have to rescale the disturbance
W_V=[[0,w_max];[0,-w_max]];
W=Polyhedron(W_V);
H_w = [0 ,1;0,-1;1,0;-1,0];
h_w = [w_max;w_max;0;0];

%Cost
Q = [1 0; 0 1e-2]; 
R = 1e-1;

% Define unitary Hyper-cube
HB_p = [1 0 ; -1 0 ; 0 1 ; 0 -1];
hB_p = 0.5*[1;1;1;1];
B_p = Polyhedron(HB_p,hB_p);

%Constraints
F = [1/4.1 0; 1/-0.2 0;0,1/5;0 1/-5; zeros(2,2)];     
G = [zeros(4,1); 1/15; 1/-15];
%Inital condition
x0=[4;0];
s0=0;
%% Offline Computation Robust Adaptive Learning MPC 
[Theta_HC0,theta_bar0,eta_0,K,P,mu,H,rho_theta0,L_B,d_bar,c,c_max]=RALMPC_Offline(Q,R,F,G,B_p,W,A_0,A_1,A_2,B_0,B_1,B_2,m,n,p);
%% Compute inital solution 
[SS_0,J_wc_0,X_bar_inital,S_inital]=get_Initalsolution(x0,s0,theta_bar0,B_p,eta_0,rho_theta0,L_B,d_bar,c,c_max,H,A_0,A_1,A_2,B_0,B_1,B_2,K,Q,R,P,F,G,m,n,p);
%% Solve RMPC
disp("Solve the RMPC")
%No Adaption to the unknown parameter
%RMPC Setting
N_RMPC=25;
M = 10;
adpation=false;
%initalize Delta
Delta=[];
%Solve the Optimization Problem
tic
[X_hat_OL_RMPC,X_bar_OL_RMPC,V_OL_RMPC,S_OL_RMPC,J_RMPC,X_RMPC,U_RMPC,time_RMPC,Theta_HC_t_RMPC]=solve_RAMPC(x0,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,P,K,F,G,d_bar,L_B,c_max,c,H,B_p,n,m,N_RMPC,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,T_s,adpation,disturbance_deter);
toc
%% Solve the RAMPC
disp("Solve the RAMPC")
%Adpation to the unkwon parameter (Paper: Linear robust adaptive model predictive control: Computational complexity and conservatism)
%RAMPC Setting
N_RAMPC = 25;
M = 10;
adpation=true;
%initalize Delta for the moving window hypercube update
for i = 1:M  
    Delta{i}=Theta_HC0;
end
%Solve the Optimization Problem
tic
[X_hat_OL_RAMPC,X_bar_OL_RAMPC,V_OL_RAMPC,S_OL_RAMPC,J_RAMPC,X_RAMPC,U_RAMPC,time_RAMPC,Theta_HC_t_RAMPC]=solve_RAMPC(x0,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,P,K,F,G,d_bar,L_B,c_max,c,H,B_p,n,m,N_RAMPC,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,T_s,adpation,disturbance_deter);
toc
%% Solve for the Optimal Solution
disp("Solve the OS")
%Setting
N_OS=100;
%Solve the Optimal Control Problem without disrubance and uncertainty (best case)
[X_OL_OS,V_OL_OS,J_OS,X_OS,U_OS,time_OS]=solve_OptimalSolution(x0,A_star,B_star,K,F,G,Q,R,n,m,N_OS,W_V);
%% Query
disp("Solve the Query")
for N=2:-2:2
    disp("N: "+N)
    %Solve the RLMPC
    %RLMPC Setting
    numberitertions=40;
    M = 0;
    adpation=false;
    %initalize Delta for the moving window hypercube update
    Delta{1}=0;
    %Set Up the inital sample set
    SS=SS_0 ;
    Q_func=J_wc_0;
    %Solve the Optimization Problem
    [X_RLMPC,U_RLMPC,S_RLMPC,J_RLMPC,X_hat_OL_RLMPC,X_bar_OL_RLMPC,V_OL_RLMPC,S_OL_RLMPC,J_wc_RLMPC,time_RLMPC,Theta_HC_t_RLMPC]=RALMPC_Convex_Online(x0,SS,Q_func,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,K,F,G,d_bar,L_B,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,T_s,numberitertions,adpation,disturbance_deter);

    %Solve the RALMPC
    %RALMPC Setting
    numberitertions=40;
    M = 10;
    adpation=true;
    %initalize Delta for the moving window hypercube update
    for i = 1:M  
        Delta{i}=Theta_HC0;
    end
    %Set Up the inital sample set
    SS=SS_0 ;
    Q_func=J_wc_0;
    %Use the convex robust safe set implementation 
    [X_RALMPC,U_RALMPC,S_RALMPC,J_RALMPC,X_hat_OL_RALMPC,X_bar_OL_RALMPC,V_OL_RALMPC,S_OL_RALMPC,J_wc_RALMPC,time_RALMPC,Theta_HC_t_RALMPC]=RALMPC_Convex_Online(x0,SS,Q_func,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,K,F,G,d_bar,L_B,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,T_s,numberitertions,adpation,disturbance_deter);
    %Plot and Save
    plotten(H,X_bar_inital,S_inital,X_RALMPC,X_RMPC,X_RAMPC,X_OS,Theta_HC0,Theta_HC_t_RALMPC,J_RMPC,J_RAMPC,J_RALMPC,J_RLMPC,J_OS,N,numberitertions,disturbance_deter)
end

