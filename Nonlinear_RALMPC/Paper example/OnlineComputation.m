clc
clear
close all
%% Load offline computed parameter
load("Parameter_Offline.mat")
p=2;
r=6;
HB_p = [1, 0; -1,0;0,1;0,-1];
hB_p = [1; 1;1;1];
B_p = Polyhedron(HB_p, hB_p);
%% Declare MPC Parameter
import casadi.*
% Prediction horizon
N=12;
% make symbolic x_bar x_hat u_bar s w  (n*(N+1)+n*(N+1)+m*N+(N+1)+2^p*N)
y=MX.sym('y',n*(N+1)+n*(N+1)+m*N+(N+1)+N);
% Number of MPC iterations
mpciterations = 100;
%% Other Optimization Options
options = optimset('Display','none',...
    'TolFun', 1e-10,...
    'MaxIter', 1000,...
    'TolConSQP', 1e-8);
%% Other Parameters
x0=[0.1;0.1];
theta_star=[1;1];
%For plotting the ellipses
[eig_U, eig_D] = eig(P);
theta_circle = linspace(0, 2*pi, 100);
circle = [cos(theta_circle); sin(theta_circle)];
%% MPC Initialization
%nonlinear constraints 
[con] = nonlinearconstraints(N, y, n, m,h_value,theta_0,theta_0,eta_0,L_B_rho,rho_theta0,L,l,c,x_s,c_xs,P,p,d_bar,delta_loc);
con_lb=[zeros(n*(N),1);zeros(n*(N),1);zeros((N),1);zeros((N),1);-inf*ones(r*N,1);-inf*ones(2,1)];
con_ub=[zeros(n*(N),1);zeros(n*(N),1);zeros((N),1);zeros((N),1);zeros(r*N,1);zeros(2,1)];
%input constraints: use lb, ub
lb=[x0;-inf*ones(n*(N),1);x0;-inf*ones(n*(N),1);-inf*ones(m*(N),1);zeros(N+1,1);zeros(N,1)];
ub=[x0;inf*ones(n*(N),1);x0;inf*ones(n*(N),1);inf*ones(m*(N),1);0;inf*ones((N),1);inf*ones((N),1)];
% objective function
obj=costfunction(N, y, x_s, u_s, Q, R,P,n,c_alpha,m);
%Option for casadi solver
opts = struct;
opts.ipopt.tol = 1e-10;
opts.ipopt.constr_viol_tol = 1e-15;
opts.ipopt.acceptable_tol = 1e-10;
opts.ipopt.max_iter =5000;
opts.ipopt.print_level = 0;  % Suppress solver output
opts.print_time = false;      % Disable timing information
%Define the problem and the solver
nlp = struct('x', y, 'f', obj, 'g', con);
solver = nlpsol('solver', 'ipopt', nlp, opts); %,'file_print_level',5

% Initial guess for input
L_Theta=eta_0*L_B_rho;
y_init=[repmat(x0,N+1,1);repmat(x0,N+1,1);zeros(m*N,1)];
for k=0:N-1
   y_init=[y_init;(1-(rho_theta0+L_Theta)^k)/(1-(rho_theta0+L_Theta))*d_bar+10^-4];
end
y_init=[y_init;0;d_bar*ones(N,1)];
% Set variables for output
x = [];
u = [];
%% Parameter Update Initialization
M=10;
theta_bar_t{1}=theta_0;
theta_hat_t{1}=theta_0;
Theta_t{1}=Theta_0;
rho_theta_t{1}=rho_theta0;
eta_t{1}=eta_0;

for j=1:M
    Delta{j}=Theta_0;
end
%% MPC Algorithm 
% initilization of measured values
xmeasure = x0;
for ii = 1:mpciterations % maximal number of iterations
    % Set initial constraint
    lb(1:n)=xmeasure;
    ub(1:n)=xmeasure;
    lb(n*(N+1)+1:n*(N+1)+n)=xmeasure;
    ub(n*(N+1)+1:n*(N+1)+n)=xmeasure;
    % Update the parameter theta_bar and theta_hat
    if ii>1
        [theta_bar_t{ii},eta_t{ii},Theta_t{ii},Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t{end},eta_t{end},Theta_0,Theta_t{end},Delta,D,options,B_p,h_value,p);
        theta_hat_t{ii}=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_t{end},theta_hat_t{end},mu,p,options,h_value);
        rho_theta_t{ii}=rho_theta0+(eta_0-eta_t{end})*L_B_rho;
    end
    % Update the constrain set and the optimization problem
    [con] = nonlinearconstraints(N, y, n, m,h_value,theta_bar_t{end},theta_hat_t{end},eta_t{end},L_B_rho,rho_theta_t{end},L,l,c,x_s,c_xs,P,p,d_bar,delta_loc);
    nlp = struct('x', y, 'f', obj, 'g', con);
    solver = nlpsol('solver', 'ipopt', nlp, opts);
    %Solve the optimization problem
    res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
    %check if the solver was succesful
    stats = solver.stats();
    if strcmp(stats.return_status, 'Solve_Succeeded')
        disp('Optimization was successful!');
    else
        disp(['Solver failed with status: ', stats.return_status,'in iteration ',ii]);
    end
    %Open loop solution
    y_OL=full(res.x); 
    X_bar_OL{ii}=reshape(y_OL(1:n*(N+1)),2,[]);
    X_hat_OL{ii}=reshape(y_OL(n*(N+1)+1:n*(N+1)+n*(N+1)),2,[]);
    U_bar_OL{ii}=y_OL(n*(N+1)+n*(N+1)+1:n*(N+1)+n*(N+1)+m*N);
    S_OL{ii}=y_OL(n*(N+1)+n*(N+1)+m*N+1:n*(N+1)+n*(N+1)+m*N+N+1);
    %Check constraints
    %[con] = nonlinearconstraints(N, y_OL, n, m,h_value,theta_bar_t{end},theta_hat_t{end},eta_t{end},L_B_rho,rho_theta_t{end},L,l,c,x_s,c_xs,P,p,d_bar,delta_loc);
    %[con, ceq] = nonlinearconstraints_paper(N,y_OL,x_s,c,c_xs,L(:,m+1:end),L(:,1:m),rho_theta_t{end},d_bar,eta_t{end},L_B_rho,theta_bar_t{end},theta_hat_t{end},delta_loc,P,n,m,B_p);
    epsilon=671.9418;
    %cost = costfunction_paper(N, y_OL, x_s, u_s, Q, R, epsilon, P,n,m)
    %cost = costfunction(N, y_OL, x_s, u_s, Q, R,P,n,c_alpha,m)
    %Close Loop
    u_cl=U_bar_OL{end}(1);
    x=[x,xmeasure];
    u=[u,u_cl];
    x_tminus=xmeasure;
    u_tminus=u_cl;
    xmeasure=simulation(x(:,end),u_cl,theta_star,D,h_value,n);
    %New initial solution for the optimization
    u_init = controller_K(u_s,x_s(1),x_s(2))*(X_bar_OL{end}(:,end)-x_s)+u_s;
    x_init = system_f(theta_bar_t{end}(1),theta_bar_t{end}(2),u_init,X_bar_OL{end}(1,end),X_bar_OL{end}(2,end));
    s_init=rho_theta_t{end}*S_OL{end}(end)+eta_t{end}*L_B_rho*S_OL{end}(end)+d_bar+max(uncertainty_w_Theta(eta_t{end},X_bar_OL{end}(1,end),X_bar_OL{end}(2,end)));
    w_init=max(uncertainty_w_Theta(eta_t{end},X_bar_OL{end}(1,end),X_bar_OL{end}(2,end)));
    y_init=[y_OL(n+1:n*(N+1));x_init;y_OL(n*(N+1)+n+1:n*(N+1)+n*(N+1));x_init;y_OL(n*(N+1)+n*(N+1)+m+1:n*(N+1)+n*(N+1)+m*N);u_init;y_OL(n*(N+1)+n*(N+1)+m*N+1+1:n*(N+1)+n*(N+1)+m*N+N+1);s_init;y_OL(n*(N+1)+n*(N+1)+m*N+N+1+1+1:end);w_init];
end
%% Plot
figure(1); hold on; axis equal;grid on;
plot([-0.1, 0.1, 0.1, -0.1, -0.1], [-0.1, -0.1, 0.1, 0.1, -0.1], 'k-', 'LineWidth', 2);
%Plot Open Loop solution
for ii=1:mpciterations
   plot(X_bar_OL{ii}(1,:),X_bar_OL{ii}(2,:))
   plot(X_bar_OL{ii}(1,:), X_bar_OL{ii}(2,:), 'ro', 'MarkerFaceColor', 'r'); % Center point
   if ii==1 || mod(ii,20)==0  
       for i=1:length(S_OL{ii})
            A_trans = eig_U * sqrt(inv(eig_D))*(S_OL{ii}(i));
            ellipse = A_trans * circle + X_bar_OL{ii}(:,i);
            plot(ellipse(1,:), ellipse(2,:), 'b', 'LineWidth', 2);
            
       end
       A_trans = eig_U * sqrt(inv(eig_D)) * (c_xs-S_OL{ii}(end));
       ellipse = A_trans * circle + x_s;
       plot(ellipse(1,:), ellipse(2,:), 'b', 'LineWidth', 3);
       plot(x_s(1), x_s(2), 'ro', 'MarkerFaceColor', 'r'); % Center point
   end
end
xlabel("x_{1,OL}")
ylabel("x_{2,OL}")
figure(2)
plot(linspace(0,mpciterations*h_value,mpciterations),x(1,:))
grid on
xlabel("time steps")
ylabel("x_{1,cl}")
figure(3)
plot(x(1,:),x(2,:))
grid on
xlabel("x_{1,cl}")
ylabel("x_{2,cl}")
figure(4)
plot(linspace(0,mpciterations*h_value,mpciterations),u)
grid on
xlabel("time steps")
ylabel("u_{cl}")
%% MPC Help functions
function [con] = nonlinearconstraints(N, y, n, m,h,theta_bar,theta_hat,eta_t,L_B_rho,rho_thetat,L,l,c,x_s,c_xs,P,p,d_bar,delta_loc)
    % Introduce the equality constraints also for the terminal state  
       x_bar=y(1:n*(N+1));
       x_hat=y(n*(N+1)+1:n*(N+1)+n*(N+1));
       u_bar=y(n*(N+1)+n*(N+1)+1:n*(N+1)+n*(N+1)+m*N);
       s=y(n*(N+1)+n*(N+1)+m*N+1:n*(N+1)+n*(N+1)+m*N+N+1);
       w=y(n*(N+1)+n*(N+1)+m*N+N+1+1:end);
       con = [];
    % dynmaics for x_bar
    for k=1:N
        x_bar_k=x_bar((k-1)*n+1:k*n);
        x_bar_new=x_bar(k*n+1:(k+1)*n);        
        u_bar_k=u_bar((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_bar_new-system_f(theta_bar(1),theta_bar(2),u_bar_k,x_bar_k(1),x_bar_k(2));
        con = [con; ceqnew];
    end
    % dynmaics for x_hat
    for k=1:N
        x_hat_k=x_hat((k-1)*n+1:k*n);
        x_hat_new=x_hat(k*n+1:(k+1)*n);
        %controller u_hat=K*(x_bar-x_hat)+u_bar
        u_bar_k=u_bar((k-1)*m+1:k*m);
        x_bar_k=x_bar((k-1)*n+1:k*n);
        u_hat_k=controller_K(u_bar_k,x_bar_k(1),x_bar_k(2))*(x_hat_k-x_bar_k)+u_bar_k;
        %dynamic constraint
        ceqnew=x_hat_new-system_f(theta_hat(1),theta_hat(2),u_hat_k,x_hat_k(1),x_hat_k(2));
        con = [con; ceqnew];
    end
    %dynamic s and upper bound w
    con_s=[];
    con_w=[];
    for k=1:N
        x_bar_k=x_bar((k-1)*n+1:k*n);
        s_k=s((k-1)+1:k);
        s_new=s(k+1:(k+1));
        w_k=w((k-1)*m+1:k*m);
        %Upper bound on the uncertainty
        w_Theta = uncertainty_w_Theta(eta_t,x_bar_k(1),x_bar_k(2));
        
        cineqnew=10^5*(max(w_Theta)^2-w_k^2);
        con_w = [con_w;  cineqnew];
        %Dynamic s
        w_delta_Theta_D=w_k+(eta_t*L_B_rho)*s_k+d_bar;
        ceqnew=-s_new+rho_thetat*s_k+w_delta_Theta_D;
        con_s = [con_s;  ceqnew];
    end
    con = [con;con_s;con_w];
    %Constraints
    for k=1:N
        u_bar_k=u_bar((k-1)*m+1:k*m);
        x_bar_k=x_bar((k-1)*n+1:k*n);
        s_k=s((k-1)+1:k);
        cineqnew=L*[u_bar_k;x_bar_k]-l+c*s_k;
        con = [con; cineqnew];
    end
    %Terminal set
    x_bar_N=x_bar(end-n+1:end);
    s_N=s(end);
    cineqnew=(x_bar_N-x_s).'*P*(x_bar_N-x_s)-(c_xs-s_N)^2;
    con = [con; cineqnew];
    cineqnew=s_N-delta_loc;
    con = [con; cineqnew];
end


function cost = costfunction(N, y, x_s, u_s, Q, R,P,n,c_alpha,m)
    % Formulate the cost function to be minimized
    cost = 0;
    x_bar=y(1:n*(N+1));
    x_hat = y(n*(N+1)+1:n*(N+1)+n*(N+1));
    u_bar=y(n*(N+1)+n*(N+1)+1:n*(N+1)+n*(N+1)+m*N);
    
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        x_hat_k=x_hat((k-1)*n+1:k*n);
        %controller u_hat=K*(x_bar-x_hat)+u_bar
        u_bar_k=u_bar((k-1)*m+1:k*m);
        x_bar_k=x_bar((k-1)*n+1:k*n);
        u_hat_k=controller_K(u_bar_k,x_bar_k(1),x_bar_k(2))*(x_hat_k-x_bar_k)+u_bar_k;
        %running kost at k
        cost = cost + runningcosts(x_hat_k, u_hat_k, x_s, u_s, Q, R);
    end
    x_bar_N=x_bar(end-n+1:end);
    cost = cost + terminalcosts( x_bar_N, x_s, P,c_alpha);
    
end

function cost = runningcosts(x, u, x_eq, u_eq, Q, R)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
end


function cost = terminalcosts(x, x_eq, P,c_alpha)
    % Introduce the terminal cost
    cost = (x-x_eq)'*P*(x-x_eq)*c_alpha;
end

function x_kplus=simulation(x_k,u_k,theta_star,D,h,n)
    w = D.V(1,:)';
    disturbance=(1-2*rand(n,1)).*w;
    x_kplus=system_f(theta_star(1),theta_star(2),u_k,x_k(1),x_k(2))+disturbance;
end

function [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,D,options,B_p,h_value,p)
    %This function updates the hypercube with Algorithm 1
    Theta_HC_tminus=Theta_HC_t;
    eta_tminus = eta_t;
    theta_bar_tminus = theta_bar_t;
    H_w=D.A;
    h_w=D.b;
    %Update Delta
    for i=1:length(Delta)-1
        Delta{i}=Delta{i+1};
    end 
    %Compute the latest Delta
    h_delta = h_w - H_w*(xmeasure - system_f(0,0,u_tminus,x_tminus(1),x_tminus(2)));
    H_delta = -H_w*system_G(x_tminus(1),x_tminus(2));
    Delta{end} = Polyhedron(H_delta,h_delta);
    Delta{end}.minHRep;
    %Compute Theta_M_t
    Theta_M_t = Theta_HC0; %Why
    for i = 1:length(Delta)
        Theta_M_t = Theta_M_t & Delta{i};
        Theta_M_t.minHRep;
    end
    H_theta_M_t=[Theta_M_t.A;Theta_HC_tminus.A];
    h_theta_M_t=[Theta_M_t.b;Theta_HC_tminus.b];
    Theta_M_t=Polyhedron(H_theta_M_t,h_theta_M_t);
    %Compute updated hypercube Theta_HC_t
        %Solve LPs
        for i = 1 : p
            e_i=[zeros(i-1,1);1;zeros(p-i,1)];
            [~,theta_min(i)]=linprog(e_i,H_theta_M_t,h_theta_M_t,[],[],[],[],options);
            [~,theta_max(i)]=linprog(-e_i,H_theta_M_t,h_theta_M_t,[],[],[],[],options);
            theta_max(i) = -theta_max(i);
        end
        eta_t = 0.5*round(max(theta_max - theta_min),5);
        %Compute the new theta_bar_t
        theta_bar_t = 0.5*(theta_max + theta_min)';
        for i=1:p   %Projection    
           if theta_bar_t(i)< theta_bar_tminus(i)-(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)-(eta_tminus-eta_t);
           elseif theta_bar_t(i)>theta_bar_tminus(i)+(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)+(eta_tminus-eta_t);             
           else
             %no projection necessary  
           end          
        end
    Theta_HC_t=theta_bar_t+eta_t*B_p;
end

function theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,mu,p,options,h_value)
%This function performs the point estimate
    theta_hat_tminus=theta_hat_t;
    %predicted state
    x_hat_1t=system_f(theta_hat_tminus(1),theta_hat_tminus(2),u_tminus,x_tminus(1),x_tminus(2));
    %Update theta_tilde
    theta_tilde_t = theta_hat_tminus + mu*system_G(x_tminus(1),x_tminus(2))*(xmeasure-x_hat_1t);
    [theta_hat_t,~] = fmincon(@(theta) norm(theta-theta_tilde_t),zeros(p,1),Theta_HC_t.A,Theta_HC_t.b,[],[],[],[],[],options);
end