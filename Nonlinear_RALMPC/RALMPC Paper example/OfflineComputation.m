clear
clc
%% Define Variables
syms x_1k x_2k x_1kplus x_2kplus u_k
syms h
syms theta_1 theta_2 eta_t 
syms v z1 z2
syms s

%% Define constants
n=2; %number of states
m=1; %number of inputs
p=2; %number of uncertain parameter
r=6; %number of constrains
h_value=0.05; %stepsize
accuracy =12; %grid accuracy
Q = 0.1*eye(n); 
R = eye(m);
rho = 0.996;
%% Define the model
x_1kplus=x_1k+h*(0.5*(1+x_1k)*u_k-x_2k*theta_1);
x_2kplus=x_2k+h*(0.5*(1-4*x_2k)*u_k+x_1k*theta_2);

% System representation
x = [x_1k; x_2k]; % State vector
theta = [theta_1; theta_2]; % Parameter vector
f = [x_1kplus; x_2kplus]; % System dynamics
f = subs(f,h,h_value);
u=u_k;                     % Input vector
% Compute Jacobian of f w.r.t. theta
G = jacobian(f, theta);
G=subs(G,h,h_value);
% Define system dimensions
m = size(u, 1);
n = size(x, 1);

% Generate MATLAB functions for system dynamics and Jacobian
matlabFunction(f, "File", "system_f");
matlabFunction(G, "File", "system_G");
% Compute the Jacobian function for A and B
compute_jacobian(f,x,u)
%% Define constraints
% State constraints
x1_min=-0.1;
x1_max=0.1;
x2_min=-0.1;
x2_max=0.1;
% Input constraints
u_min=-2;
u_max=2;
% Constraint matrices
L_u = kron(eye(m), [-1; 1]);
l_u = [-u_min; u_max];
L_x = kron(eye(n), [-1; 1]);
l_x = [-x1_min; x1_max; -x2_min; x2_max];
L_x=[-10,0;10,0;0,-10;0,10];
l_x=ones(2^n,1);
L_u=[-0.5;0.5];
l_u=ones(2^m,1);
% Combined state and input constraints
L = [L_u, zeros(size(L_u, 1), size(L_x, 2));
     zeros(size(L_x, 1), size(L_u, 2)), L_x];
l = [l_u; l_x];

% Compute upper and lower bound constraints
con_u_lb = u_min;
con_u_ub = u_max;
con_x_lb = [x1_min;x2_min];
con_x_ub = [x1_max;x2_max];
%% Parametric Uncertainty
HB_p = [1, 0; -1,0;0,1;0,-1];
hB_p = [1; 1;1;1];
B_p = Polyhedron(HB_p, hB_p);
theta_0 = [1.01;0.99];
eta_0 = 0.01;
Theta_0 = theta_0 + eta_0 * B_p;
%% Additive Disturbance 
D =0.5*10^-4*B_p;
%% Set up optimization problem
fprintf("Computing P and K\n");
% Define decision variables
X = sdpvar(n);
Y_0 = sdpvar(m, n);
Y_1 = sdpvar(m, n);
Y_2 = sdpvar(m, n);
Y_3 = sdpvar(m, n);
Y_4 = sdpvar(m, n);
Y_5 = sdpvar(m, n);
Y_6 = sdpvar(m, n);
Y_7 = sdpvar(m, n);
Y_8 = sdpvar(m, n);
Y_9 = sdpvar(m, n);
delta=[v,z1,z2,v^2,z1^2,z2^2,v*z1,v*z2,z1*z2];
delta_fun=@(v,z1,z2) [v,z1,z2,v^2,z1^2,z2^2,v*z1,v*z2,z1*z2];

% Define optimization settings
ops = sdpsettings('solver', 'sdpt3', 'verbose', 1);        
% Define objective function
obj = -log(det(X));

% Compute constraints for the Linear Matrix Inequalities (LMI)
fprintf("\tCompute constraints for the LMI\n");
[con]=construct_lmi_controller(rho,con_u_lb,con_u_ub,con_x_lb,con_x_ub,L,Theta_0,accuracy,X,Y_0,Y_1,Y_2,Y_3,Y_4,Y_5,Y_6,Y_7,Y_8,Y_9,delta_fun,m);

% Solve the optimization problem
fprintf("\t Solve the SDP\n");
diagnostics = optimize(con, obj,ops);
%check the feasibility
X = value(X);
Y = value(Y_0)+delta(1)*value(Y_1)+delta(2)*value(Y_2)+delta(3)*value(Y_3)+delta(4)*value(Y_4)+delta(5)*value(Y_5)+delta(6)*value(Y_6)+delta(7)*value(Y_7)+delta(8)*value(Y_8)+delta(9)*value(Y_9);
fprintf("\t Check the solution\n");
con=check_lmi_constrain(X,Y,con_u_lb,con_u_ub,con_x_lb,con_x_ub,Theta_0,3*accuracy,rho,L,m,z1,z2,v);
if ~any(con==0) && min(real(eig(X)))>0
    fprintf('SDP successfully solved!\n');
else
    fprintf('Solving SDP failed. Please try again.\n');
end

% Obtain solution for P and K
P = inv(X);
K = Y / X;
matlabFunction(K,"File","controller_K");

%% Compute rho and delta_loc
fprintf("Searching delta_loc\n");
%delta_loc=22.81;
delta_loc=sqrt(23.8689);
accuracy=10;
[rho_theta0,delta_loc,Psi]=compute_contraction(delta_loc,P,L_x,l_x,L_u,l_u,con_u_lb,con_u_ub,con_x_lb,con_x_ub,accuracy,theta_0);
%% Compute L_B_rho
 fprintf("Computing L_B_rho\n");
%Compute the set Psi: Using results from the previous section
% Define the range and number of grid points for each state
%Set up the SDP
L_B_rho=compute_L_B(Psi,P,B_p,n);

%% Compute d_bar
fprintf("Solving SDP for d_bar");
d_bar=computing_dbar(P,D);
%% Compute c_j
fprintf("Solving SDP for c_j\n");
gamma=sdpvar(1);
c=[];
con=[];
fprintf("\tIterate over constraints\n");
Psi_subset = Psi(:, randperm(size(Psi, 2), round(size(Psi, 2)*0.5)));
for i=1:size(L,1)
    obj=gamma;
    c_j=compute_c_j(Psi,P,L,l,n,i   );
    c=[c;c_j];
end
%% Compute the parameter update gain mu
fprintf("Compute the parameter update gain mu\n");
options = optimset('Display','off');
[test,fval] = fmincon(@(x)  1/getDsq(x),[1;1],...
    L_x,l_x,[],[],[],[],[],options);
mu = floor(fval); 
%% Terminal set
fprintf("Compute the terminal set\n");
x_s=[0;0]   ;
u_s=0;
fprintf("\tCheck if the Point is a steady state \n");
if all(x_s== system_f(theta_0(1),theta_0(2),u_s,x_s(1),x_s(2)))
    disp("It is a steady state!")
end
%c_xs
c_xs=min([-(L*[u_s;x_s]-l)./c;delta_loc]);
%Check the scalar condition 
fprintf("\tCheck the scalar condition \n")
L_theta0=eta_0*L_B_rho;
w_theta0=[];
for i=1:size(Theta_0.V,1)
    w_theta0(i)=sqrt((system_G(x_s(1),x_s(2))*Theta_0.V(i,:).').'*P*(system_G(x_s(1),x_s(2))*Theta_0.V(i,:).'));
end
w_theta0D_s=eta_0*max(w_theta0)+d_bar;

if rho_theta0+L_theta0+c_xs*w_theta0D_s<=1
    disp("The terminal set condition is satisfied")
else
    disp("The terminal set condition is not satisfied")
end
%Compute the factor c_alpha
K_s= controller_K(u_s,x_s(1),x_s(2));
c_alpha=max(eig(Q+K_s.'*R*K_s))/((1-(rho_theta0+L_theta0)^2)*min(eig(P)));
%% Compute the function for the uncertainty discription 
w_Theta =[];
for j=1:size(B_p.V,1)
    L_Theta=eta_t*L_B_rho;
    G_z=system_G(z1,z2);
    w_Theta=[w_Theta;eta_t*sqrt((G_z*B_p.V(j,:).').'*P*(G_z*B_p.V(j,:).'))];
end
matlabFunction(w_Theta, "File", "uncertainty_w_Theta");
%% Save the parameter
Y_0=value(Y_0);Y_1=value(Y_1);Y_2=value(Y_2);Y_3=value(Y_3);Y_4=value(Y_4);Y_5=value(Y_5);Y_6=value(Y_6);Y_7=value(Y_7);Y_8=value(Y_8);Y_9=value(Y_9);
save("Parameter_Offline.mat","P","n","m","p","r","h_value","Q","R","rho_theta0","delta_loc","mu","c","c_xs","c_alpha","Theta_0","eta_0","theta_0","L","l","x_s","u_s","D","d_bar","L_B_rho","B_p","Y_0","Y_1","Y_2","Y_3","Y_4","Y_5","Y_6","Y_7","Y_8","Y_9")
%% Help functions
function compute_jacobian(f,x,u)
A=jacobian(f,x);
A=subs(A);
matlabFunction(A,"File","jacobian_A");
B=jacobian(f,u);
B=subs(B);
matlabFunction(B,"File","jacobian_B");
end

function [con]=construct_lmi_controller(rho,con_u_lb,con_u_ub,con_x_lb,con_x_ub,L,Theta_0,accuracy,X,Y_0,Y_1,Y_2,Y_3,Y_4,Y_5,Y_6,Y_7,Y_8,Y_9,delta,m)
    con=[];     
    for v=linspace(con_u_lb,con_u_ub,accuracy)
        for z1=linspace(con_x_lb(1),con_x_ub(1),accuracy) 
            for z2=linspace(con_x_lb(2),con_x_ub(2),accuracy)
                del=delta(v,z1,z2);
                Y=Y_0+del(1)*Y_1+del(2)*Y_2+del(3)*Y_3+del(4)*Y_4+del(5)*Y_5+del(6)*Y_6+del(7)*Y_7+del(8)*Y_8+del(9)*Y_9;
                for j=1:size(Theta_0.V(:,1),1)
                    theta_1=Theta_0.V(j,1);
                    theta_2=Theta_0.V(j,2);
                    A = jacobian_A(theta_1,theta_2,v);
                    B = jacobian_B(z1,z2);
                    con=[con;[rho^2*X,(A*X+B*Y)';(A*X+B*Y),X]>=0];
                end      
            end
        end
    end
    for v=linspace(con_u_lb,con_u_ub,2)
        for z1=linspace(con_x_lb(1),con_x_ub(1),2) 
            for z2=linspace(con_x_lb(2),con_x_ub(2),2)
                for j=1:size(L,1)
                    del=delta(v,z1,z2);
                    Y=Y_0+del(1)*Y_1+del(2)*Y_2+del(3)*Y_3+del(4)*Y_4+del(5)*Y_5+del(6)*Y_6+del(7)*Y_7+del(8)*Y_8+del(9)*Y_9;
                    con=[con;[1,(L(j,m+1:end)*X+L(j,1:m)*Y);(L(j,m+1:end)*X+L(j,1:m)*Y).',X]>=0];
                end
            end
        end
    end
end

function con=check_lmi_constrain(X,Y,con_u_lb,con_u_ub,con_x_lb,con_x_ub,Theta_0,accuracy,rho,L,m,z1,z2,v)
    con=[];     
    accuracy=accuracy;
    for v_value=linspace(con_u_lb,con_u_ub,accuracy)
        for z1_value=linspace(con_x_lb(1),con_x_ub(1),accuracy) 
            for z2_value=linspace(con_x_lb(2),con_x_ub(2),accuracy)
                Y_value=double(subs(Y,[z1,z2,v],[z1_value,z2_value,v_value]));
                for j=1:size(Theta_0.V(:,1),1)
                    theta_1=Theta_0.V(j,1);
                    theta_2=Theta_0.V(j,2);
                    A = jacobian_A(theta_1,theta_2,v_value);
                    B = jacobian_B(z1_value,z2_value);
                    con=[con;min(real(eig([rho^2*X,(A*X+B*Y_value).';(A*X+B*Y_value),X])))>=0]; 
                end
                
            end
        end
    end
    for v_value=linspace(con_u_lb,con_u_ub,accuracy)
        for z1_value=linspace(con_x_lb(1),con_x_ub(1),accuracy) 
            for z2_value=linspace(con_x_lb(2),con_x_ub(2),accuracy)
                Y_value=double(subs(Y,[z1,z2,v],[z1_value,z2_value,v_value]));
                for j=1:size(L,1)
                    con=[con;min(real(eig([1,(L(j,m+1:end)*X+L(j,1:m)*Y_value);(L(j,m+1:end)*X+L(j,1:m)*Y_value).',X])))>=0];
                end
            end
        end
    end
end

function [rho_theta0,delta_loc,Psi]=compute_contraction(delta_loc,P,L_x,l_x,L_u,l_u,con_u_lb,con_u_ub,con_x_lb,con_x_ub,numPoints,theta_0)
    rho_theta0=0;
    num_points=500000;
    Psi=[];
    % Precompute grid points
    x1 = linspace(con_x_lb(1), con_x_ub(1), numPoints);
    x2 = linspace(con_x_lb(2), con_x_ub(2), numPoints);
    v = linspace(con_u_lb, con_u_ub, numPoints);
    
    % Create a 3D grid of states
    [X1, X2] = ndgrid(x1, x2);
    [V] = ndgrid(v);
    
    % Reshape into vectors
    states = [X1(:), X2(:)];
    states_z = states(1:1:end, :);
    inputs = [V(:)];
    
    % Precompute sizes
    numStates_z = size(states_z, 1);
    numInputs = size(inputs, 1);
    
    
    % Iterate over subset of states
    for i = 1:numStates_z
        z = states_z(i, :).';  % Single z point
        % Randomly sample x in ||x-z||_P<=delta_loc
        x_samples = sampleEllipsoid(P, z, delta_loc, num_points,L_x,l_x);
        x_samples = [x_samples,z];
        for j=1:numInputs 
            v=inputs(j);
            z_plus = system_f(theta_0(1),theta_0(2),v,z(1),z(2)); %compute the next state z_plus

            for k=1:length(x_samples) %iteration of the samples in x_sample
                x=x_samples(:,k);
                u = controller_K(v,z(1),z(2)) * (x - z)+v; %compute the controller u

                if all(L_u * u <= l_u, 1) %Check input constraints
                   x_plus = system_f(theta_0(1),theta_0(2),u,x(1),x(2));%compute the next state x_plus

                   if all(L_x*z_plus-l_x<=0) && all(L_x*x_plus-l_x<=0) %Check if the next state are in the constrains 
                            Psi=[Psi,[x;z;v]]; %Add a value to set Psi
                        if ((x - z).' * P * (x - z)) ~= 0 
                            %Compute the contraction rate
                            rho_theta0=max([rho_theta0,sqrt((x_plus - z_plus).' * P * (x_plus - z_plus))/sqrt((x - z).' * P * (x - z))]);
                        end  
                   end
                end
            end
        end

    end
end

function L_B_rho=compute_L_B(Psi,P,B_p,n)
con=[];
    for i=1:size(Psi,2)
            z=Psi(n+1:2*n,i);
            x=Psi(1:n,i);
            G_x = system_G(x(1),x(2));
            G_z = system_G(z(1),z(2));
            for j=1:size(B_p.V,1)
                if (x-z)~=0
                    con=[con;(((G_x-G_z)*B_p.V(j,:).').'*P*((G_x-G_z)*B_p.V(j,:).'))/((x-z).'*P*(x-z))];
                end
            end
    end
    L_B_rho=max(con);
end

function con=computing_dbar(P,D)
con=[];
    for i=1:size(D.V,1)
        con=[con;D.V(i,:)*P*D.V(i,:).'];
    end
    con=max(con);
end

function c_j=compute_c_j(Psi,P,L,l,n,j)
con=[];
    for i=1:size(Psi,2)
            z=Psi(n+1:2*n,i);
            x=Psi(1:n,i);
            v=Psi(2*n+1:end,i);
            u=controller_K(v,z(1),z(2))*(x-z)+v;
            h_x=L(j,:)*[u;x]-l(j,:);
            h_z=L(j,:)*[v;z]-l(j,:);
            V_delta=sqrt((x-z).'*P*(x-z));
            con=[con;(h_x-h_z)/V_delta];
    end
    c_j=max(con);
end

function [y_sq] = getDsq(X)
    y = system_G(X(1),X(2));
    y_sq = norm(y)^2; 
end

function x = sampleEllipsoid(P, z, delta_loc, num_points,L_x,l_x)
    % Generate random points inside an ellipsoid ||x - z||_P <= delta_loc
    % P: Shape matrix (positive definite)
    % z: Center of the ellipsoid (column vector)
    % delta_loc: Scaling factor
    % num_points: Number of points to sample
    % Output: points (each column is a sampled point)

    n = length(z);  % Dimension of space
    
    %Eigenvalues and transformation matrix 
    [Q,D] = eig(P);
    
    %generate random number
    randnumber=rand(n,num_points);
    %scalings for cosinus and sinus
    a=randnumber(1,:)*delta_loc/sqrt(D(1,1));
    b=randnumber(1,:)*delta_loc/sqrt(D(2,2));

    x= Q * [a .* cos(2*pi*randnumber(2,:)); b .* sin(2*pi*randnumber(2,:))] + z;

    %Check constraints
    x=x(:,all(L_x*x<=l_x));
end