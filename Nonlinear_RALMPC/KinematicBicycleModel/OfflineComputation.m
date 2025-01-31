clear 
clc
%% Load Track data
track=load("L_track_barc.mat");
%% Define Variables
syms s_k s_kplus e_yk e_ykplus e_psik e_psikplus v_k v_kplus delta_k delta_kplus %states
syms a_k u_deltak %inputs
syms h l_f l_r kappa theta
%% Define the model
beta_k=atan((l_r/(l_r+l_f))*tan(delta_k));
s_kplus=s_k+h*v_k*(cos(e_psik+beta_k)/(1-e_yk*kappa));
e_ykplus=e_yk+h*v_k*sin(e_psik+beta_k);
e_psikplus=e_psik+h*(v_k/l_f*sin(beta_k)-kappa*v_k*(cos(e_psik+beta_k)/(1-e_yk*kappa)));
v_kplus=v_k+h*theta*a_k;
delta_kplus=delta_k+h*u_deltak;
%system
x=[s_k;e_yk;e_psik;v_k;delta_k];
u=[u_deltak;a_k];
f=[s_kplus;e_ykplus;e_psikplus;v_kplus;delta_kplus];
%size
m=size(u,1);
n=size(x,1);
%generate a matlab function
matlabFunction(f,"File","system_f");
%% Define constants
h_value=0.01;
l_f_value=0.125;%[m]
l_r_value=0.125;%[m]

accuracy=2;

Q=eye(n);
R=eye(m);
rho=0.9995;
%% Define constraints
%State
s_min=0;
s_max=sum(track.cl_segs(:,1));
e_ymin=-track.track_width/2;
e_ymax=track.track_width/2;
e_psimin=-25/180*pi;
e_psimax=25/180*pi;
v_min=0.5;
v_max=3.5;
delta_min=-0.4;
delta_max=0.4;
%Input
u_deltamin=-3;
u_deltamax=3;
a_min=-1.3;
a_max=3;
%Matrix
L_u = kron(eye(m),[-1;1]);
l_u=[-u_deltamin;u_deltamax;
    -a_min;a_max];
L_x= kron(eye(n),[-1;1]);
l_x=[-s_min;s_max;
    -e_ymin;e_ymax;
    -e_psimin;e_psimax;
    -v_min;v_max;
    -delta_min;delta_max];
% Combined state and input constraints
L = [L_u,zeros(size(L_u,1),size(L_x,2));
   zeros(size(L_x,1),size(L_u,2)),L_x];
l = [l_u;l_x];
%upper and lower bound constraints
con_u_lb=-l_u(1:2:end);
con_u_ub=l_u(2:2:end);
con_x_lb=-l_x(1:2:end);
con_x_ub=l_x(2:2:end);
%% Parametric Uncertainty
HB_p = [1 ; -1 ];
hB_p = [1;1];
B_p = Polyhedron(HB_p,hB_p);
theta_0=1;
eta_0=0.01;
Theta_0=theta_0+eta_0*B_p;
%% Additive Disturbance 
A_d=kron(eye(n),[-1;1]);
b_d=0.001*[1;1;
    e_ymax;e_ymax;
    e_psimax;e_psimax;
    v_max;v_max;
    delta_max;delta_max];
D = Polyhedron(A_d,b_d);
%% Compute the Jacobian function for A and B
compute_jacobian(f,x,u,h,l_f,l_r,h_value,l_f_value,l_r_value)
%% Set up optimization problem
% Define decision variables
X = sdpvar(n);
Y = sdpvar(m,n);
% Define objective
obj = -log(det(X));
%Gridding
[A_grid,B_grid]=gridding(track,con_x_lb,con_x_ub,Theta_0,accuracy);
%Compute constraints for the LMI
con=construct_lmi_constrain(X,Y,A_grid,B_grid,rho,L,n,m);
%Optimization
diagnostics = optimize(con,obj);
%% Obtain solution for P and K
X = value(X);
Y = value(Y);
P = inv(X);
K = Y/X;
check=check_lmi_constrain(X,Y,A_grid,B_grid,rho,L,n,m);
%% Compute rho and delta_loc
% Define the range and number of grid points for each state
numPoints = 10; % Number of points per state
delta_loc=10000;
x1 = linspace(s_min, s_max, numPoints);
x2 = linspace(e_ymin, e_ymax, numPoints);
x3 = linspace(e_psimin, e_psimax, numPoints);
x4 = linspace(v_min, v_max, numPoints);
x5 = linspace(delta_min, delta_max, numPoints);
u1 = linspace(u_deltamin, u_deltamax, numPoints);
u2 = linspace(a_min, a_max, numPoints);
% Create a 3D grid of states
[X1, X2, X3,X4,X5] = ndgrid(x1, x2, x3,x4,x5);
[U1,U2]=ndgrid(u1,u2);
% Reshape grid matrices into state vectors
states = [X1(:), X2(:), X3(:), X4(:), X5(:)];
inputs = [U1(:), U2(:)];
numStates = size(states, 1);
numInputs = size(inputs,1);
% Compute differences for all possible (x, z) pairs
cond=true;
while cond
    for i = 1:numStates
        for j = 1:numStates
            for k=1:numInputs
                x = states(i, 1:n).'; % First point
                z = states(j, 1:n).'; % Second point
                v = inputs(k, :).';
                u=K*(x-z)+v;
                % Compute the difference x - z
                diff = x - z;
                norm_diff_P = sqrt(diff.' * P * diff);            
                if norm_diff_P<=delta_loc && ~any(L_u*u-l_u>0)
                    kappa=eval_kappa(z(1),track);
                    z_plus=system_f(v(2),z(5),z(3),z(2),h_value,kappa,l_f_value,l_r_value,z(1),theta_0,v(1),z(4));      
                    x_plus=system_f(u(2),x(5),x(3),x(2),h_value,kappa,l_f_value,l_r_value,x(1),theta_0,u(1),x(4));
                    norm_diff_P_plus=sqrt((x_plus-z_plus).' * P * (x_plus-z_plus));
                    rho*norm_diff_P
                    pause(0.5)
                    if norm_diff_P_plus>rho*norm_diff_P 
                        cond=false;
                        break
                    end    
                end
            end
            if cond==false
                break
            end
        end
        if cond==false
            break
        end
    end
    if cond==false
       cond=true;
       delta_loc=delta_loc-100;
    else
        break
    end
end
%% Compute d_bar
gamma=sdpvar(1);
obj=gamma;
con=construct_lmi_d(gamma,P,D);
diagnostics = optimize(con,obj);
d_bar=value(gamma);
%% Help function 
function [A_grid,B_grid]=gridding(track,con_x_lb,con_x_ub,Theta_0,accuracy)
idx=[find(track.cl_segs(:,2)==0,1);find(track.cl_segs(:,2)~=0)];
distance=[0;cumsum(track.cl_segs(:,1))];
count=0;
for k=1:length(Theta_0.V)
    for i=1:length(idx)
        s=mean(distance(idx(i):idx(i)+1));
        kappa=eval_kappa(s,track);
        for e_y=linspace(con_x_lb(2),con_x_ub(2),accuracy)
            for e_psi=linspace(con_x_lb(3),con_x_ub(3),accuracy)
                for v=linspace(con_x_lb(4),con_x_ub(4),2) %enters affine
                    for delta=linspace(con_x_lb(5),con_x_ub(5),accuracy)
                        A = jacobian_A(delta,e_psi,e_y,kappa,v);
                        B=jacobian_B(Theta_0.V(k));
                        count=count+1;
                        A_grid(:,:,count)=A;
                        B_grid(:,:,count)=B;
                    end
                end
            end    
        end
    end
end    
end

function compute_jacobian(f,x,u,h,l_f,l_r,h_value,l_f_value,l_r_value)
A=jacobian(f,x);
A=subs(A,[h;l_f;l_r],[h_value;l_f_value;l_r_value]);
matlabFunction(A,"File","jacobian_A");
B=jacobian(f,u);
B=subs(B,[h;l_f;l_r],[h_value;l_f_value;l_r_value]);
matlabFunction(B,"File","jacobian_B");
end

function kappa=eval_kappa(s,track)
index=find(cumsum(track.cl_segs(:,1))>s,1);
r=track.cl_segs(index,2);
if s==sum(track.cl_segs(:,1))
    r=track.cl_segs(end,2);
end
if r==0
    kappa=0;
else
    kappa=1/r;
end

end

function check=check_lmi_constrain(X,Y,A_grid,B_grid,rho,L,n,m)
check=true;
for i=1:size(A_grid,3)
    if any(round(det([rho^2*X,(A_grid(:,:,i)*X+B_grid(:,:,i)*Y).';(A_grid(:,:,i)*X+B_grid(:,:,i)*Y),X]),20)<0)
        check=false;
    end  
end
for j=1:size(L,1)
    if any(det([1,(L(j,m+1:end)*X+L(j,1:m)*Y);(L(j,m+1:end)*X+L(j,1:m)*Y).',X])<0)
        check=false;
    end
end
end
function con=construct_lmi_constrain(X,Y,A_grid,B_grid,rho,L,n,m)
con=[];
for i=1:size(A_grid,3)
    con=[con;[rho^2*X,(A_grid(:,:,i)*X+B_grid(:,:,i)*Y).';(A_grid(:,:,i)*X+B_grid(:,:,i)*Y),X]>=1e-10];
end
for j=1:size(L,1)
    con=[con;[1,(L(j,m+1:end)*X+L(j,1:m)*Y);(L(j,m+1:end)*X+L(j,1:m)*Y).',X]>=1e-10];
end
end

function con=construct_lmi_d(gamma,P,D)
con=[];
    for i=1:size(D.V,1)
        con=[con;gamma^2-D.V(i,:)*P*D.V(i,:).'>=0];
    end
    con=[con;gamma>=0];
end