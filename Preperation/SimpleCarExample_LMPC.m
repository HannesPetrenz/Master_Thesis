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
%Weightings
Q=diag([1;50]);
%Constrain bounds
X=[0,inf;0,inf];
U=[-2,2];
%% Compute the first feasible solution
[x_feasible,u_feasible]=intialfeasible(100,n,m,x_init,A,B);
x_feasible=reshape(x_feasible,2,[]);
figure
plot(x_feasible(1,:),'-o')
hold on
%% Initialize Safe Set and Q-funtion
x_cl{1} = x_feasible;                          % Safe set vector: Vector collecting the state of the performed iterations
u_cl{1} = u_feasible;                         % Safe Input set vector: Vector collecting the input of the performed iterations[ Qfun] = ComputeCost(x,Q,xF);
Qfun=ComputeCost(x_feasible,Q,xF);
IterationCost{1} = Qfun;
%% Now run the LMPC
% Pick number of iterations to run
Iterations = 20;
% Run the LMPC
[ x_LMPC, u_LMPC, x_cl, u_cl, IterationCost, SS] = LMPC(x_init, x_cl, u_cl,xF, IterationCost, A, B, Q, N, Iterations, X, U,n,m);
%% Compute Optimal Solution
[ x_opt, u_opt, cost_opt] = solve_CFTOCP(x_init,xF, 300, Q, A, B, X, U,n,m,x_cl{end},u_cl{end});
%% Plotting
plot(x_opt(1,:),'LineWidth',2)
hold on
grid on
for i=1:Iterations
    plot(x_cl{i}(1,:))
end
hold off 

h = legend('First Feasible Solution','Optimal Solution','interpreter','latex');
xlabel('$$x_1$$','interpreter','latex','fontsize',20);
ylabel('$$x_2$$','interpreter','latex','fontsize',20);
set(h,'fontsize',15);
figure(2)
plot(cellfun(@(v)v(1),IterationCost))
hold on
yline(cost_opt(1))
grid on 
xlabel('$$Iteration step$$','interpreter','latex','fontsize',20);
ylabel('$$Cost$$','interpreter','latex','fontsize',20);
xlim([1,Iterations])
%% Generate initial feasible trajectory 
function [x0,u0]=intialfeasible(N,n,m,x_init,A,B)
    x0=zeros(n*(N+1),1);
    u0=zeros(m*(N),1);
    x0(1:n)=x_init;
    for i=1:N
        if i==1
            u_=1.5;    
        elseif i==67
            u_=-1;
        elseif i==69
            u_=-0.5;
        else
            u_=0;
        end
        x0(n*i+1:n*i+2)=A*x0(n*i-1:n*i)+B*u_;
        u0((i-1)*m+1:i*m)=u_;
    end
end
%% Compute Cost 
function [ Qfun] = ComputeCost(x,Q,xF)
    Cost=zeros(1,size(x,2));
    for i = (size(x,2))-1:-1:1
        Cost(i) =Cost(i+1)+runningcosts(x(:,i),xF,Q);    
    end
    Qfun=Cost;
end
function cost=runningcosts(x_k,xF,Q)
    cost=(x_k-xF)'*Q*(x_k-xF);
end
function cost = costfunction(y,N,n,xF,Q)
    % Formulate the cost function to be minimized
    
    cost = 0;
    x = y(1:n*(N+1));
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N+1
        x_k=x(n*(k-1)+1:n*k);
        cost = cost + runningcosts(x_k,xF,Q);
    end
end
%% LMPC Problem
function [ x_LMPC, u_LMPC, x_cl_out, u_cl_out, IterationCost_out, SS] = LMPC(x0, x_cl, u_cl,xF, IterationCost, A, B, Q, N, Iterations, X, U,n,m)
j = 1;
SS = x_cl{1};
uSS = u_cl{1};
Qfun= IterationCost{1};
IterationCost_out{1} = IterationCost{1};
x_cl_out{1} = x_cl{1};
u_cl_out{1} = u_cl{1};
while (j <= Iterations)
    SSQfun = Polyhedron([SS', Qfun']);
    SSQfun.computeHRep()
    SSQfun.computeVRep()
    SS = SSQfun.V(:,1:2)';
    Qfun = SSQfun.V(:, end)';

    t = 1;        % Initialize time
    x_LMPC = x0;  % Initial condition

    tollerance = 10^(-3);
    exitFlag = 0;
    while ( exitFlag == 0 )
        if t==0
            x_init=x_cl_out{j};
            u_init=u_cl_out{j};
        else
            [minValue,closestIndex] = min(vecnorm(x_cl_out{j}-repmat(x_LMPC(:,end),[1 length(x_cl_out{j})])));
            x_init=x_cl_out{j}(closestIndex:end,:);
            u_init=u_cl_out{j}(closestIndex:end,:);
        end
        [ x_OL, uPred ] = FTOCP(x_LMPC(:,t),xF,x_cl_out{j},u_cl_out{j}, N,n,m, Q, Qfun, SS, A, B, X, U);
        u_LMPC(:,t) = uPred;
        x_LMPC(:,t+1) = dynamic(x_LMPC(:,t),u_LMPC(:,t),A,B);

        t = t + 1;
        if (x_LMPC(:,t)-xF)'*(x_LMPC(:,t)-xF) <(tollerance)
            exitFlag = 1;
        end
    end
    % Now save the data, update cost and safe set.
    x_cl_out{j+1} = x_LMPC;
    u_cl_out{j+1} = u_LMPC.';
    IterationCost_out{j+1} = ComputeCost(x_LMPC,Q,xF);
    
    SS   = [SS, x_LMPC];
    Qfun = [Qfun, IterationCost_out{j+1}];

    % increase Iteration index and restart
    j = j + 1;
    if j <= Iterations
        clear x_LMPC
        clear u_LMPC
    end
    
end
end

function [ x_OL, uPred ] = FTOCP(x_t,xF,x_init,u_init, N,n,m, Q, Qfun, SS, A, B, X, U)
import casadi.*
%make symbolicm, optimization variable
y=MX.sym('y',(N)*m+(N+1)*n+length(Qfun));
%input constraints: use lb, ub
lb=[repmat(X(:,1),(N+1),1);repmat(U(:,1),(N),1);zeros(length(Qfun),1)];    %[X,U,Constraint the multipliers to be positive]
ub=[repmat(X(:,2),(N+1),1);repmat(U(:,2),(N),1);inf*ones(length(Qfun),1)]; %[X,U,Constraint the multipliers to be positive]
%initial state constraint: use lb, ub
lb(1:n)=x_t;
ub(1:n)=x_t;
%objective function
obj=costfunction_FTOCP(y,N,n,m,xF,Q,Qfun);
%nonlinear constraints (dynamics,Terminal point in the convex hull,Must be convex combination --> sum to 1)
con_bound=zeros((N)*n+n+1,1);
con_lb=con_bound;
con_ub=con_bound;
con=nonlinearconstraints_FTOCP(N, y, n, m,A,B,Qfun,SS) ;
%construct the nonlinear programming problem
nlp = struct('x', y, 'f', obj, 'g', con);
options = struct;
options.ipopt.print_level=0;
options.ipopt.tol=1e-10;
solver = casadi.nlpsol('solver','ipopt',nlp,options);
% Set initial guess and initial constraint
if length(x_init)<=N+1
    y_init=[reshape(x_init,1,[]).';repmat(x_init(:,end),N-length(x_init)+1,1);u_init;repmat(u_init(end,:),N-length(x_init)+1,1);1;zeros(length(Qfun)-1,1)];
else
    x_init=reshape(x_init,1,[]).';
    y_init=[x_init(1:n*(N+1),1);u_init(1:m*N,1);1;zeros(length(Qfun)-1,1)];
end
res = solver('x0' , y_init,... % solution guess
                 'lbx', lb,...           % lower bound on x
                 'ubx', ub,...           % upper bound on x
                 'lbg', con_lb,...           % lower bound on g
                 'ubg', con_ub);             % upper bound on g
if solver.stats.return_status~="Solve_Succeeded"
    fprintf('No solution is feasible \n')
end
y_OL=full(res.x); 
x_OL=reshape(y_OL(1:n*(N+1)),n,[]);
u_OL=reshape(y_OL(n*(N+1)+1:n*(N+1)+m*N),m,[]);
lambda_OL=y_OL(n*(N+1)+m*N+1:end);
uPred=u_OL(1);
end

function cost = costfunction_FTOCP(y,N,n,m,xF,Q,Qfun)
    % Formulate the cost function to be minimized
    lambda=y((N)*m+(N+1)*n+1:end);
    cost = 0;
    x = y(1:n*(N+1));
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N+1
        x_k=x(n*(k-1)+1:n*k);
        cost = cost + runningcosts(x_k,xF,Q);
    end
    %terminal cost
    cost=cost+Qfun*lambda;
end
function [con] = nonlinearconstraints_FTOCP(N, y, n, m,A,B,Qfun,SS) 

   % Introduce the nonlinear constraints also for the terminal state   
   x = y(1:n*(N+1));
   x_s = y(n*N+1:n*(N+1));%x_N
   lambda=y(n*(N+1)+m*N+1:end);
   u = y(n*(N+1)+1:n*(N+1)+m*N);
   con = [];

   % system dynamics constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);       
        u_k=u((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - dynamic(x_k,u_k,A,B);
        con = [con; ceqnew];
    end
   % Terminal Constraint: enforce predicted state in the convex safe set
   ceqnew=x_s-SS*lambda; % Terminal point in the convex hull
   con = [con; ceqnew];
   ceqnew=1-ones(1,length(Qfun))*lambda; % Must be convex combination --> sum to 1
   con = [con; ceqnew];
end
%% Optimal Control Problem
function [ x_cl, u_cl,cost] = solve_CFTOCP( x_t,xF, N, Q, A, B, X, U,n,m,x_init,u_init)
[lb,ub,con_lb,con_ub,solver]=CreateFTOCP(X,U,N,m,n,x_t,xF,Q,A,B);
% Set initial guess and initial constraint
if length(x_init)<=N+1
    y_init=[reshape(x_init,1,[]).';repmat(x_init(:,end),N-length(x_init)+1,1);u_init;repmat(u_init(end,:),N-length(x_init)+1,1)];
else
    x_init=reshape(x_init,1,[]).';
    y_init=[x_init(1:n*(N+1),1);u_init(1:m*N,1)];
end
lb(1:n)=x_t;
ub(1:n)=x_t;
%Solving the open loop problem
res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
% check whether the problem is feasible
if solver.stats.return_status~="Solve_Succeeded"
    fprintf('No solution is feasible \n')
end
y_OL=full(res.x); 
x_cl=reshape(y_OL(1:n*(N+1)),n,[]);
u_cl=reshape(y_OL(n*(N+1)+1:end),m,[]);
cost= ComputeCost(x_cl,Q,xF);
end
function [lb,ub,con_lb,con_ub,solver]=CreateFTOCP(X,U,N,m,n,x_init,xF,Q,A,B)
import casadi.*
%make symbolicm, optimization variable
y=MX.sym('y',(N)*m+(N+1)*n);
%input constraints: use lb, ub
lb=[repmat(X(:,1),(N+1),1);repmat(U(:,1),(N),1)];
ub=[repmat(X(:,2),(N+1),1);repmat(U(:,2),(N),1)];
%initial state constraint: use lb, ub
lb(1:n)=x_init;
ub(1:n)=x_init;
%objective function
obj=costfunction(y,N,n,xF,Q);
%nonlinear constraints (dynamics)
con_bound=zeros((N)*n,1);
con_lb=con_bound;
con_ub=con_bound;
con=nonlinearconstraints(N, y, n, m,A,B);
%construct the nonlinear programming problem
nlp = struct('x', y, 'f', obj, 'g', con);
options = struct;
options.ipopt.print_level=0;
options.ipopt.tol=1e-10;
solver = casadi.nlpsol('solver','ipopt',nlp,options);
end
function [con] = nonlinearconstraints(N, y, n, m,A,B) 

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
   ceqnew=x_s - dynamic(x((N-1)*n+1:N*n), u(end-m+1:end), A,B);%without terminal constraint
   con = [con; ceqnew];
end
%% Dynamics
function x=dynamic(x0,u,A,B)
x=A*x0+B*u;
end