    clear 
clc
close all
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
F = [1/1.1 0; 1/-0.1 0;0,1/5;0 1/-5; zeros(2,2)];     
G = [zeros(4,1); 1/5; 1/-5];
%Inital condition
x0=[1;0];
s0=0;
%% Offline Computation Robust Adaptive Learning MPC 
[Theta_HC0,theta_bar0,eta_0,K,P,mu,H,rho_theta0,L_B,d_bar,c,c_max]=RALMPC_Offline(Q,R,F,G,B_p,W,A_0,A_1,A_2,B_0,B_1,B_2,m,n,p);
%% Compute inital solution 
[SS_0,J_wc_0]=get_Initalsolution(x0,s0,theta_bar0,B_p,eta_0,rho_theta0,L_B,d_bar,c,c_max,H,A_0,A_1,A_2,B_0,B_1,B_2,K,Q,R,P,F,G,m,n,p);
%% Solve the RAMPC (without Learning the terminal set and cost),
%RAMPC Setting
N = 14;
M = 10;
%initalize Delta for the moving window hypercube update
for i = 1:M  
    Delta{i}=Theta_HC0;
end
tic
[X_hat_OL,X_bar_OL,V_OL,S_OL,J_wc,X,U,time,eta_t,rho_theta_t,theta_bar_t,theta_hat_t,Theta_HC_t]=solve_RAMPC(x0,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,P,K,F,G,d_bar,L_B,c_max,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,T_s);
toc
%% Plot
figure(2)
for t=1:length(X_hat_OL)
    for k=1:length(X_hat_OL{t})
        x_bar_OL=X_bar_OL{t};
        s_OL=S_OL{t};
        X_t=Polyhedron(H,s_OL(k)*ones(size(H,1),1)+H*x_bar_OL(:,k));
        plot(X_t)
        hold on
    end
    plot(x_bar_OL(1,:),x_bar_OL(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20)
    grid on
    plot(X(1),X(2),"Marker",'.', 'MarkerSize', 20)
end
xlabel("x_{1}")
ylabel("x_{2}")
title("X_{k|t}")
% figure(3)
% if t==0
%     Theta_HC0.plot('color', 'red');
% else
%     hold on
%     Theta_HC_t.plot('color', 'lightblue');
% end
