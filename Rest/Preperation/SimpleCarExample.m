close all
clear all
clc
import casadi.*
%% Declaration
% System dynmics
A=eye(2,2)+[0 1;0 0];
B=[0;1];
% System dimensions
n = size(A,1); % state dimension
m = size(B,2); % input dimension
% Initial condition
x_init=[0;0];
% Terminal State
xF=[100;0];
% Horizon 
N =3;
% compute inital guess
[x0,u0]=intialfeasible(N,n,m,x_init,A,B);
%Creating the optimazation problem
[lb,ub,con_lb,con_ub,solver]=genOptProblem(N,m,n,x_init,xF,A,B);
%% MPC scheme 
% Set variables for output
xmeasure = x_init;
x = xmeasure;
u = [];
it=0;

while norm(x(:,end)-xF)>=0.01 && it<=200
    % Set initial guess and initial constraint
    y_init=[x0;u0];
    lb(1:n)=xmeasure;
    ub(1:n)=xmeasure;
    %Solving the open loop problem
    res = solver('x0' , y_init,... % solution guess
                 'lbx', lb,...           % lower bound on x
                 'ubx', ub,...           % upper bound on x
                 'lbg', con_lb,...           % lower bound on g
                 'ubg', con_ub);             % upper bound on g
    % check whether the problem is feasible
    if solver.stats.return_status~="Solve_Succeeded"
        fprintf('No solution is feasible \n')
        break
    end
    % optimal solution
    y_OL=full(res.x); 
    x_OL=y_OL(1:n*(N+1));
    u_OL=y_OL(n*(N+1)+1:end);
    % Update closed-loop system (apply first control move to system)
    xmeasure = x_OL(n+1:2*n);
    % Store closed loop data
    x = [ x, xmeasure ];
    u = [ u, u_OL(1:m) ];
    % Compute initial guess for next time step
    u0 = [u_OL(m+1:end);zeros(m,1)];
    x0 = [x_OL(n+1:end);zeros(n,1)];
    % Update iteration
    it=it+1;
end
%% Compute Optimal solution
N =200;
% compute inital guess
[x0,u0]=intialfeasible(N,n,m,x_init,A,B);
%Creating the optimazation problem
[lb,ub,con_lb,con_ub,solver]=genOptProblem(N,m,n,x_init,xF,A,B);
%Solve the optimal problem
xmeasure = x_init;
% Set initial guess and initial constraint
y_init=[x0;u0];
lb(1:n)=xmeasure;
ub(1:n)=xmeasure;
res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
%optimal solution
y_opt=full(res.x); 
x_opt=reshape(y_opt(1:n*(N+1)),2, []);
u_opt=y_opt(n*(N+1)+1:end);
%% Plot
figure(1)
subplot(3,1,1)
plot(x(1,:))
hold on
plot(x_opt(1,:),'*')
hold off
grid on
legend("optimal solution","initial solution")
subplot(3,1,2)
plot(x(2,:))
hold on
plot(x_opt(2,:),'*')
hold off
grid on
legend("optimal solution","initial solution")
subplot(3,1,3)
plot(u)
hold on
plot(u_opt)
grid on
legend("optimal input","initial input")
%% Dynamics
function x=dynamic(x0,u,A,B)
x=A*x0+B*u;
end
%% Object function
function cost = costfunction(y,N,n,xF)
    % Formulate the cost function to be minimized
    
    cost = 0;
    x = y(1:n*(N+1));
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N+1
        x_k=x(n*(k-1)+1:n*k);
        cost = cost + runningcosts(x_k,xF,n);
    end
end

function cost=runningcosts(x_k,xF,n)
    cost=(x_k-xF)'*diag([1,50])*(x_k-xF);%if_else(x_k==xF, 1, 0);
end


%% Constraints
function [con] = nonlinearconstraints(N, y, n, m,A,B,xF) 

   % Introduce the nonlinear constraints also for the terminal state   
   x = y(1:n*N);
   x_s = y(n*N+1:n*(N+1));%x_N
   u = y(n*(N+1)+1:end);
   con = [];

   % system dynamics constraints along prediction horizon
    for k=1:N-1
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k,u_k,A,B);
        con = [con; ceqnew];
    end
   % Terminal constraint
   ceqnew=x_s - dynamic(x((N-1)*n+1:N*n), u(end-m+1:end), A,B);%xF - dynamic(x((N-1)*n+1:N*n), u(end-m+1:end), A,B);
   con = [con; ceqnew];
end
%% Generate initial feasible trajectory 
function [x0,u0]=intialfeasible(N,n,m,x_init,A,B)
    x0=zeros(n*(N+1),1);
    u0=zeros(m*(N),1);
    x0(1:n)=x_init;
    for i=1:N
        if i==1
            u_=1;    
        elseif i==6
            u_=-1;
        else
            u_=0;
        end
        x0(n*i+1:n*i+2)=A*x0(n*i-1:n*i)+B*u_;
        u0((i-1)*m+1:i*m)=u_;
    end
end
%% Generate Optimitation Problem
function [lb,ub,con_lb,con_ub,solver]=genOptProblem(N,m,n,x_init,xF,A,B)
import casadi.*
%make symbolicm, optimization variable
y=MX.sym('y',(N)*m+(N+1)*n);
%input constraints: use lb, ub
lb=[zeros(n*(N+1),1);-2*ones(m*(N),1)];
ub=[inf*ones(n*(N+1),1);2*ones(m*(N),1)];

%initial state constraint: use lb, ub
lb(1:n)=x_init;
ub(1:n)=x_init;
%terminal constraint:
%objective function
obj=costfunction(y,N,n,xF);
%nonlinear constraints (dynamics, terminal constraint, artificial steady-state constraint)
con_bound=zeros((N)*n,1);
con_lb=con_bound;
con_ub=con_bound;
con=nonlinearconstraints(N, y, n, m,A,B,xF);
%construct the nonlinear programming problem
nlp = struct('x', y, 'f', obj, 'g', con);
options = struct;
options.ipopt.print_level=0;
options.ipopt.tol=1e-10;
solver = casadi.nlpsol('solver','ipopt',nlp,options);
end