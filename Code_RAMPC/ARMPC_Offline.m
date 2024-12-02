% OFFLINE COMPUTATIONS
clear all
close all
clc
options = optimset('Display','off');
%% Constraints on state and input
x_1min = -0.1;
x_1max = 1.1;
x_2min = -5;
x_2max = 5;
u_min = -5;
u_max = 5;

% F and G from the constraints above
F = [1/x_1max 0; 1/x_1min 0;0,1/x_2max;0 1/x_2min; zeros(2,2)];     
G = [zeros(4,1); 1/u_max; 1/u_min];
fgk = size(F,1);  % Total number of constraints we have (6 here)
%% Parameters and uncertainties
p = 2;      % parameter dimension
n = 2;      % state dimension
m = 1;      % input dimension

% sampling time for first order discretization
T_s = 0.1;
% uncertainties reduction to satisfy the condition on terminal region existance
unc_reduction = 1;

mass = 1;
spring = 1;
damper = 0.2;
par_uncertainty = 0.5;
disturbances = 0.2;  

% What we consider as parameters
% p1 = k
% p2 = c
% d' = d/m --> D changes

p1_unc = par_uncertainty*(spring/mass);
p2_unc = par_uncertainty*(damper/mass);
w_max = T_s/m*disturbances; 
%disturbance set 

W_V=[[0,w_max];[0,-w_max]];
W=Polyhedron(W_V);
H_w = [0 ,1;0,-1;1,0;-1,0];
h_w = [w_max;w_max;0;0];
%parameter set
H_p = [1 0; -1 0; 0 1; 0 -1];
h_p = unc_reduction*[1;1;1;1];
Theta = Polyhedron(H_p,h_p);
% unitary Hyper-cube
HB_p = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
hB_p = 0.5*[1;1;1;1;1;1];
B_p = Polyhedron(HB_p,hB_p);
%% System matrices + Euler discretization
A_0 = [1 T_s; -T_s*spring/mass 1-T_s*damper/mass];
A_1 = p1_unc*[0 0; T_s 0];
A_2 = p2_unc*[0 0; 0 T_s];
B_0 = [0; T_s/mass];
B_1 = zeros(n,m);
B_2 = zeros(n,m);
%% Costs
Q = [1 0; 0 1e-2]; 
R = 1e-1;
%% Start time
t_start = tic;
%% Solve for P and K (Assumption 2)
rho_PK = 0.75; %set contraction rate \in[0,1) (hyperparameter), see [https://github.com/Gedlex/nonlinear-robust-MPC/tree/main] for automatic choice
[P,K] = get_PK_LMIs(Q,R,Theta,rho_PK,A_0,A_1,A_2,B_0,B_1,B_2);
%% A_K matrices at vertices of $\Theta$
for i = 1:2^p
    vert_Theta = Theta.V(i,:)';
    eval(['A_K' num2str(i) '=getA_K_full(vert_Theta,K,A_0,A_1,A_2,B_0,B_1,B_2);']);
end 
%% Compute set X_0 as rho contractive. Compute the set P
rho = rho_PK; 
% symmetric constraints to build X_0
x_1max_symm = -x_1min;
% F from the constraint above
F_symm = [1/x_1max_symm 0; 1/x_1min 0;F(3:end,:)];
%also shift input constraint to make it symmetric
G_symm = [zeros(4,1); 1/(u_max-1); 1/u_min];
q=size(F_symm,1);
disp('Computing X_0 as rho contractive set....');
disp('Iteration step:');
X_0 = lambda_contractive(rho,F_symm,G_symm,A_K1,A_K2,A_K3,A_K4,K);

X_0.minHRep;
X_0.minVRep;
Hx = X_0.A;
hx = X_0.b;
c_0 = length(Hx(:,1));
% normalize constraints 
for i = 1:c_0
    Hx(i,:) = Hx(i,:)/hx(i);
    hx(i) = 1;
end
% Plot X_0's verteces and compute max x2 in X_0 just to plot it decently later
X0_2max = -Inf;
%% Plot contractive set used for tbe
figure(1);
plot(X_0)
vert_X_0=con2vert(X_0.A,X_0.b);%converts to vertex representation 
hold on
for l = 1:length(Hx(:,1))
    x0 = vert_X_0(l,:)';
    plot(x0(1),x0(2),'og','LineWidth',1)
    if x0(2) > X0_2max
        X0_2max = x0(2);
    end
end

%% Compute constants for robust MPC
%rho_0
array_rho = [];
[A_rho,B_rho] = getAB(zeros(p,1),A_0,A_1,A_2,B_0,B_1,B_2);
A_K = A_rho + B_rho*K;
for i = 1:c_0
    [~,fval] = linprog(-Hx(i,:)*A_K, Hx, ones(size(Hx,1),1),[],[],[],[],options);
    array_rho = [array_rho; -f  val];
end
rho_0 = max(array_rho);
% L_theta_0
array_L_theta = [];
for j = 1:(2^p)
    [A_Lj,B_Lj] = getABunc(Theta.V(j,:),A_1,A_2,B_1,B_2);
    A_K = A_Lj + B_Lj*K;
    for i = 1:c_0
        [~,fval] = linprog(-Hx(i,:)*A_K, Hx, ones(size(Hx,1),1),[],[],[],[],options);
        array_L_theta = [array_L_theta; -fval];
    end
end
L_Theta_0 = max(array_L_theta);
% c_j and c_max
c_j = [];
init = ones(1,c_0);
lb_cj = zeros(c_0,1);
ub_cj = +inf*ones(c_0,1);
fun = @(x) sum(x(1:c_0));
for j = 1:fgk
    con_cj = @(x) function_con_cj(x,F+G*K,Hx,j,c_0);
    [arg,~] = fmincon(fun, init,[],[],[],[],lb_cj,ub_cj,con_cj,options);
    c_j = [c_j; max(arg)];
end
cmax = max(c_j);
%alternative computation with LP
%c_j_2=[];
%for j=1:fgk
%    [~,temp]=linprog(-F(j,:)-G(j,:)*K,Hx,ones(size(Hx,1),1));
%   c_j_2=[c_j_2;-temp]; 
%end

% dbar
array_dbar = [];
for i = 1:(n)
    dbar = Hx*W.V(i,:)';
    array_dbar = [array_dbar; dbar];
end
dbar = max(array_dbar);
%% check terminal set condition
%condition is checked here also for the steady state [1;0]
cond = rho_0+L_Theta_0+cmax*dbar;
ratio = (1-rho_0)/(L_Theta_0+cmax*dbar);
eta=2;%hard coded here recall: theta\in[-1,1]; and Theta_HC=[-0.5,0.5]^2
L_B=L_Theta_0/eta;
x_ss = [1;0];
%added computation to take care of terminal for x_s\neq 0 
u_ss = x_ss(1)/mass*spring*[1-par_uncertainty,1+par_uncertainty];% two extreme inputs
%possible maginitudes of model mismatch at steady-state
w_max_s=[];
for l=1:2
w_max_s=[w_max_s; Hx*((A_1+A_2)*x_ss+(B_1+B_2)*u_ss(l));...
        Hx*((A_1-A_2)*x_ss+(B_1-B_2)*u_ss(l));...
        -Hx*((A_1+A_2)*x_ss+(B_1+B_2)*u_ss(l));...
        -Hx*((A_1-A_2)*x_ss+(B_1-B_2)*u_ss(l));...
        ];
end
% 
f_theta=[];
for j=1:size(F,1)
    for l=1:2%exteme input
     f_theta=[f_theta;(1-F(j,:)*x_ss-G(j,:)*u_ss(l))/c_j(j)];%shifted by steady-state (Prop. 4)
    end
end
f_min=min(f_theta); 
%
%cond should be <1 for terminal set conditin to hold
cond=(rho_0+eta*L_B+max(w_max_s)+dbar)/f_min
ratio = f_min*(1-rho_0)/(eta*L_B+dbar+max(w_max_s))
%% compute mu for RLS
[test,fval] = fmincon(@(x) 1/getDsq(x,A_1,A_2,B_1,B_2),[1;1;1],...
    [F G],ones(fgk,1),[],[],[],[],[],options);
mu_gain = floor(fval);
%% measure time and display things

disp(['Condition: ',num2str(cond)]);
if cond > 1
    disp(['Uncertainties we can handle only this much uncertainty: ',num2str(unc_reduction*ratio*100),'%']);
    disp(' ');
elseif cond <= 1
    disp(['Uncertainties: ',num2str(unc_reduction*100),'%']);
    disp(' ');
end
disp(['Number of sides of X_0_LPV: ',num2str(c_0)]);
t_offline = toc(t_start);
disp(['Time for offline computations: ',num2str(t_offline),' s']); 
%% Save offline things 
save('ARMPC_Offline');
%%
%% FUNCTIONS
function [A_K] = getA_K_full(p,K,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2 ;
    B = B_0 + p(1)*B_1 + p(2)*B_2 ;
    A_K = A + B*K;
end

function [A,B] = getAB(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function [A,B] = getABunc(p,A_1,A_2,B_1,B_2)
    A = p(1)*A_1 + p(2)*A_2;
    B = p(1)*B_1 + p(2)*B_2;
end
 
function [y_sq] = getDsq(X,A_1,A_2,B_1,B_2)
    n = length(A_1(1,:));   % state dimension
    m = length(B_1(1,:));   % input dimension
    y = [A_1*X(1:n) + B_1*X(n+1:n+m), A_2*X(1:n) + B_2*X(n+1:n+m)];
    y_sq = norm(y)^2; 
end

function [c,ceq] = function_con_cj(x,L,Hx,j,c_0)
    % Extract vector h
    temp = x(1:c_0);
    h = reshape(temp,[1,c_0]);
    c = [];
    ceq = h*Hx - L(j,:);
end
 