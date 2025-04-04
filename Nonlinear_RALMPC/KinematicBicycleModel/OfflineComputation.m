clear
clc
%% Load Track data
track=load("L_track_barc.mat");

%% Define Variables
syms s_k s_kplus e_yk e_ykplus e_psik e_psikplus v_k v_kplus delta_k delta_kplus % States
syms a_k u_deltak % Inputs
syms h l_f l_r kappa theta

%% Define the model
% Compute the slip angle beta_k
beta_k = atan((l_r/(l_r + l_f)) * tan(delta_k));

% State update equations
s_kplus = s_k + h * v_k * (cos(e_psik + beta_k) / (1 - e_yk * kappa));
e_ykplus = e_yk + h * v_k * sin(e_psik + beta_k);
e_psikplus = e_psik + h * (v_k / l_f * sin(beta_k) - kappa * v_k * (cos(e_psik + beta_k) / (1 - e_yk * kappa)));
v_kplus = v_k + h * theta * a_k;
delta_kplus = delta_k + h * u_deltak;

% System representation
x = [s_k; e_yk; e_psik; v_k; delta_k]; % State vector
u = [u_deltak; a_k]; % Input vector
f = [s_kplus; e_ykplus; e_psikplus; v_kplus; delta_kplus]; % System dynamics

% Compute Jacobian of f w.r.t. theta
G = jacobian(f, theta);

% Define system dimensions
m = size(u, 1);
n = size(x, 1);

% Generate MATLAB functions for system dynamics and Jacobian
matlabFunction(f, "File", "system_f");
matlabFunction(G, "File", "system_G");

%% Define constants
h_value = 0.01;
l_f_value = 0.125; % [m]
l_r_value = 0.125; % [m]
accuracy = 10;
Q = eye(n);
R = eye(m);
rho = 0.9995;

%% Define constraints
% State constraints
s_min = 0;
s_max = sum(track.cl_segs(:,1));
e_ymin = -track.track_width / 2;
e_ymax = track.track_width / 2;
e_psimin = -25 / 180 * pi;
e_psimax = 25 / 180 * pi;
v_min = 0.5;
v_max = 3.5;
delta_min = -0.4;
delta_max = 0.4;

% Input constraints
u_deltamin = -3;
u_deltamax = 3;
a_min = -1.3;
a_max = 3;

% Constraint matrices
L_u = kron(eye(m), [-1; 1]);
l_u = [-u_deltamin; u_deltamax; -a_min; a_max];
L_x = kron(eye(n), [-1; 1]);
l_x = [-s_min; s_max; -e_ymin; e_ymax; -e_psimin; e_psimax; -v_min; v_max; -delta_min; delta_max];

% Combined state and input constraints
L = [L_u, zeros(size(L_u, 1), size(L_x, 2));
     zeros(size(L_x, 1), size(L_u, 2)), L_x];
l = [l_u; l_x];

% Compute upper and lower bound constraints
con_u_lb = -l_u(1:2:end);
con_u_ub = l_u(2:2:end);
con_x_lb = -l_x(1:2:end);
con_x_ub = l_x(2:2:end);

%% Parametric Uncertainty
HB_p = [1; -1];
hB_p = [1; 1];
B_p = Polyhedron(HB_p, hB_p);
theta_0 = 1;
eta_0 = 0.01;
Theta_0 = theta_0 + eta_0 * B_p;
%% Additive Disturbance 
A_d=kron(eye(n),[-1;1]);
b_d=0.0001*[1;1;
    e_ymax;e_ymax;
    e_psimax;e_psimax;
    v_max;v_max;
    delta_max;delta_max];
D = Polyhedron(A_d,b_d);
%% Compute the Jacobian function for A and B
compute_jacobian(f, x, u, h, l_f, l_r, h_value, l_f_value, l_r_value);

%% Set up optimization problem
fprintf("Computing P and K\n");
% Define decision variables
X = sdpvar(n);
Y = sdpvar(m, n);
% Define optimization settings
ops = sdpsettings('verbose',0);
% Define objective function
obj = -log(det(X));

% Compute system matrices via gridding
fprintf("\tCompute system matrices via gridding\n");
[A_grid, B_grid] = gridding(track, con_x_lb, con_x_ub, Theta_0, accuracy);

% Compute constraints for the Linear Matrix Inequalities (LMI)
fprintf("\tCompute constraints for the LMI\n");
con = construct_lmi_constrain(X, Y, A_grid, B_grid, rho, L, n, m);

% Solve the optimization problem
fprintf("\t Solve the SDP\n");
diagnostics = optimize(con, obj,ops);
if strcmp(diagnostics.info,'Successfully solved (SeDuMi)')
    fprintf('SDP successfully solved!\n');
else
    fprintf('Solving SDP failed. Please try again.\n');
end
%% Obtain solution for P and K
X = value(X);
Y = value(Y);
P = inv(X);
K = Y / X;

% Check the validity of LMI constraints
check = check_lmi_constrain(X, Y, A_grid, B_grid, rho, L, n, m);
%% Compute rho and delta_loc
fprintf("Gridding to find delta_loc\n");
% Define the range and number of grid points for each state
numPoints = accuracy; % Number of points per state
delta_loc = 1000;
cond = true;

fprintf("\tPrecompute grid points\n");
% Precompute grid points
idx = [find(track.cl_segs(:,2) == 0, 1); find(track.cl_segs(:,2) ~= 0)];
distance = [0; cumsum(track.cl_segs(:,1))];
x1 = mean([distance(idx), distance(idx+1)].');
x2 = linspace(e_ymin, e_ymax, numPoints);
x3 = linspace(e_psimin, e_psimax, numPoints);
x4 = linspace(v_min, v_max, numPoints);
x5 = linspace(delta_min, delta_max, numPoints);
v1 = linspace(u_deltamin, u_deltamax, numPoints/2);
v2 = linspace(a_min, a_max, numPoints/2);

% Create a 3D grid of states
[X1, X2, X3, X4, X5] = ndgrid(x1, x2, x3, x4, x5);
[V1, V2] = ndgrid(v1, v2);

% Reshape into vectors
states = [X1(:), X2(:), X3(:), X4(:), X5(:)];
states_z = states(1:1:end, :);
inputs = [V1(:), V2(:)];

% Precompute sizes
numStates = size(states, 1);
numStates_z = size(states_z, 1);
numInputs = size(inputs, 1);

% Precompute Cholesky decomposition
L_chol = chol(P);

% Transform all states once
states_tilde = L_chol * states.';

% Build KD-Tree once (outside loop)
stateTree = KDTreeSearcher(states_tilde.');

fprintf("\tIterate over the Grid Points\n");
% Iterate over theta values
for k = 1:size(Theta_0.V, 1)
    theta = Theta_0.V(k);

    % Iterate over subset of states
    for i = 1:numStates_z
        z = states_z(i, :).';  % Single z point
        z_tilde = L_chol * z;  % Transform z

        % Find nearby states x using KD-Tree
        idx = rangesearch(stateTree, z_tilde.', delta_loc);
        x_tilde = states_tilde(:, idx{1}); 
        x = states(idx{1}, :).';  % Extract corresponding x

        % Compute u for all valid (x, z, v)
        u_xz = K * (x - z);  % Matrix operation

        % Expand `inputs` across `x` to efficiently compute `U`
        V=kron(inputs.', ones(1, size(u_xz, 2)));
        U = repmat(u_xz, 1, numInputs) + V;

        % Check constraint L_u * u - l_u <= 0
        valid_mask = all(L_u * U <= l_u, 1);

        % Extract valid indices
        valid_idx = find(valid_mask);

        % Ensure valid_x and valid_U have the same dimensions
        valid_U = U(:, valid_idx);
        valid_x = repmat(x, 1, numInputs);  % Repeat `x` for each input combination
        valid_x = valid_x(:, valid_idx);  % Select only valid entries
        valid_v =V(:,valid_idx);
        % Iterate over valid (u, x)
        for j = 1:size(valid_U, 2)
            u_ = valid_U(:, j);
            x_ = valid_x(:, j);
            v_ = valid_v(:, j);
            % Compute next states
            z_plus = system_f(v_(2), z(5), z(3), z(2), h_value, eval_kappa(z(1), track), l_f_value, l_r_value, z(1), theta, v_(1), z(4));
            x_plus = system_f(u_(2), x_(5), x_(3), x_(2), h_value, eval_kappa(x_(1), track), l_f_value, l_r_value, x_(1), theta, u_(1), x_(4));

            % Check contraction condition
            if sqrt((x_plus - z_plus).' * P * (x_plus - z_plus)) - rho * sqrt((x_ - z).' * P * (x_ - z)) > 0
                disp("Condition not satisfied. Decrease the size of delta_loc");
                cond = false;
                break
            end    
        end
    end
    if cond==false
        break
    end
end
rho_theta0=rho;
%% Compute L_B_rho
fprintf("Solving SDP for L_B_rho\n");
%Compute the set Psi: Using results from the previous section
 % Iterate over subset of states
fprintf("\tDefining the set Psi\n");
for i = 1:numStates_z
    z = states_z(i, :).';  % Single z point
    z_tilde = L_chol * z;  % Transform z
    % Find nearby states x using KD-Tree
    idx = rangesearch(stateTree, z_tilde.', delta_loc);
    x_tilde = states_tilde(:, idx{1}); 
    x = states(idx{1}, :).';  % Extract corresponding x
    % Compute u for all valid (x, z, v)
    u_xz = K * (x - z);  % Matrix operation

    % Expand `inputs` across `x` to efficiently compute `U`
    V=kron(inputs.', ones(1, size(u_xz, 2)));
    U = repmat(u_xz, 1, numInputs) + V;
    % Check constraint L_u * u - l_u <= 0
    valid_mask = all(L_u * U <= l_u, 1);

    % Extract valid indices
    valid_idx = find(valid_mask);
    
    % Ensure valid_x and valid_U have the same dimensions
    valid_u = U(:, valid_idx);
    valid_x = repmat(x, 1, numInputs);  % Repeat `x` for each input combination
    valid_x = valid_x(:, valid_idx);  % Select only valid entries
    valid_v =V(:,valid_idx);
    valid_z=repmat(z,1,size(valid_v,2));
    Psi{i}=[valid_z;valid_v;valid_x;valid_u];
end
%Set up the SDP
gamma=sdpvar(1);
obj=gamma;
fprintf("\tConstruct the LMI\n");
con=construct_lmi_L_B(gamma,Psi,P,Theta_0,h_value,n,m);
fprintf("\tSolve the SDP\n");
diagnostics = optimize(con,obj,ops);
L_B_rho=value(gamma);
if strcmp(diagnostics.info,'Successfully solved (SeDuMi)')
    fprintf('SDP successfully solved!\n');
else
    fprintf('Solving SDP failed. Please try again.\n');
end
%% Compute d_bar
fprintf("Solving SDP for d_bar");
gamma=sdpvar(1);
obj=gamma;
fprintf("\tConstruct the LMI\n");
con=construct_lmi_d(gamma,P,D);
fprintf("\tSolve the SDP\n");
diagnostics = optimize(con,obj,ops);
d_bar=value(gamma);
if strcmp(diagnostics.info,'Successfully solved (SeDuMi)')
    fprintf('SDP successfully solved!\n');
else
    fprintf('Solving SDP failed. Please try again.\n');
end
%% Compute c_j
fprintf("Solving SDP for c_j\n");
gamma=sdpvar(1);
c=[];
con=[];
fprintf("\tIterate over constraints\n");
for i=1:size(L,1)
    obj=gamma;
    con=[(L(i,m+1:end)+L(i,1:m)*K).'*(L(i,m+1:end)+L(i,1:m)*K)-gamma*P<=0;gamma>=0];
    fprintf("\tSolving LMI number %d \n",i);
    diagnostics = optimize(con,obj,ops);
    if strcmp(diagnostics.info,'Successfully solved (SeDuMi)')
        fprintf('\tSDP successfully solved!\n');
    else
        fprintf('\tSolving SDP failed. Please try again.\n');
    end
    c=[c;value(gamma)];
end
%% Compute the parameter update gain mu
fprintf("Compute the parameter update gain mu\n");
options = optimset('Display','off');
[test,fval] = fmincon(@(x) 1/getDsq(x,h_value),[1;1],...
    L_u,l_u,[],[],[],[],[],options);
mu = floor(fval);
%% Terminal set
fprintf("Compute the terminal set\n");
s_s=track.cl_segs(1,1)+track.cl_segs(2,1)+track.cl_segs(3,1)/2;
x_s=[s_s;0;0;1;0];
u_s=[0;0];
fprintf("\tCheck if the Point is a steady state \n");
if all(x_s== system_f(u_s(2), x_s(5), x_s(3), x_s(2), h_value, eval_kappa(x_s(1), track), l_f_value, l_r_value, x_s(1), theta, u_s(1), x_s(4)))
    disp("It is a steady state!")
end
%c_xs
c_xs=min([-(L*[u_s;x_s]-l)./c;delta_loc]);
%Check the scalar condition 
fprintf("\tCheck the scalar condition \n")
L_theta0=eta_0*L_B_rho;
w_theta0=[];
for i=1:size(Theta_0.V,2)
    w_theta0(i)=sqrt((system_G(u_s(2),h_value)*Theta_0.V(i)).'*P*(system_G(u_s(2),h_value)*Theta_0.V(i)));
end
w_theta0D_s=eta_0*max(w_theta0)+d_bar;

if rho_theta0+L_theta0+c_xs*w_theta0D_s<=1
    disp("The terminal set condition is satisfied")
else
    disp("The terminal set condition is not satisfied")
end

%% Save the results to a mat file
fprintf("\tSave the results\n")
save("RALMPC_OfflineCompuation",'K','P','rho_theta0','L_B_rho','delta_loc','mu','d_bar',"c_xs","L","l",'theta_0','eta_0','Theta_0')













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

function con=construct_lmi_L_B(gamma,Psi,P,Theta,h_value,n,m)
con=[];
    for i=1:length(Psi)
        for k=1:size(Psi{i},2)
            z=Psi{i}(1:n,k);
            v=Psi{i}((n+1):(n+m),k);
            x=Psi{i}((n+m+1):(n+m+n),k);
            u=Psi{i}((n+m+n+1):end,k);
            G_x = system_G(u(2),h_value);
            G_v = system_G(v(2),h_value);
            for j=1:size(Theta.V,1)
                con=[con;[0>=gamma^2*(x-z).'*P*(x-z)-((G_x-G_v)*Theta.V(j)).'*P*((G_x-G_v)*Theta.V(j))]];
            end
        end
    end
    con=[con;gamma>=0];
end

function [y_sq] = getDsq(X,h)
    y = system_G(X(2),h);
    y_sq = norm(y)^2; 
end