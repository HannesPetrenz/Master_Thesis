function [X_bar_OL,U_bar_OL,J,X,U]=solve_OptimalSolution(x0,x_s, u_s, Q, R,theta_star,L,l,r,N,mpciterations,m,n)
import casadi.*
%Option for casadi solver
opts = struct;
opts.ipopt.tol = 1e-10;
opts.ipopt.constr_viol_tol = 1e-15;
opts.ipopt.acceptable_tol = 1e-10;
opts.ipopt.max_iter =5000;
opts.ipopt.print_level = 0;  % Suppress solver output
opts.print_time = false;      % Disable timing information
% make symbolic x_bar x_hat u_bar s w  (n*(N+1)+n*(N+1)+m*N+(N+1)+2^p*N)
y=MX.sym('y',n*(N+1)+m*N);
%nonlinear constraints 
[con] = nonlinearconstraints(x_s,N, y, n, m,theta_star,L,l);
con_lb=[zeros(n*(N),1);-inf*ones(r*N,1);zeros(n,1)];
con_ub=[zeros(n*(N),1);zeros(r*N,1);zeros(n,1)];
%input constraints: use lb, ub
lb=[x0;-inf*ones(n*(N),1);-inf*ones(m*(N),1)];
ub=[x0;inf*ones(n*(N),1);inf*ones(m*(N),1)];
% objective function
obj= costfunction(N, y, x_s, u_s, Q, R,n,m);
%construct the optimization problem 
nlp = struct('x', y, 'f', obj, 'g', con);
solver = nlpsol('solver', 'ipopt', nlp, opts);
%compute inital solution
y_init=zeros(size(y,2));
%initialize more variables
t=0;
xmeasure = x0;
J{1}=0;
X=[];
U=[];
while t<mpciterations
    % Set initial constraint
    lb(1:n)=xmeasure;
    ub(1:n)=xmeasure;
    %Solve the optimization problem
    res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
    %check if the solver was succesful
    stats = solver.stats();
    if strcmp(stats.return_status, 'Solve_Succeeded')
        %disp('Optimization was successful!');
    else
        disp(['Solver failed with status: ', stats.return_status]);
    end
    %Open loop solution
    y_OL=full(res.x); 
    X_bar_OL{t+1}=reshape(y_OL(1:n*(N+1)),2,[]);
    U_bar_OL{t+1}=y_OL(n*(N+1)+1:n*(N+1)+m*N);
    %Close Loop
    u_cl=U_bar_OL{end}(1);
    X=[X,xmeasure];
    U=[U,u_cl];
    x_tminus=xmeasure;
    u_tminus=u_cl;
    xmeasure=simulation(X(:,end),u_cl,theta_star);
    %New initial solution for the optimization
    u_init = controller_K(u_s,x_s(1),x_s(2))*(X_bar_OL{end}(:,end)-x_s)+u_s;
    x_init = system_f(theta_star(1),theta_star(2),u_init,X_bar_OL{end}(1,end),X_bar_OL{end}(2,end));
    y_init=[y_OL(n+1:n*(N+1));x_init;y_OL(n*(N+1)+m+1:n*(N+1)+m*N);u_init];
    %Compute cost function
    J{end+1}=J{end}+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
    %Next step
    t=t+1;
end
end

%% Help functions
function [con] = nonlinearconstraints(x_s,N, y, n, m,theta_star,L,l)
    % Introduce the equality constraints also for the terminal state  
       x_bar=y(1:n*(N+1));
       u_bar=y(n*(N+1)+1:n*(N+1)+m*N);
       con = [];
    % dynmaics for x_bar
    for k=1:N
        x_bar_k=x_bar((k-1)*n+1:k*n);
        x_bar_new=x_bar(k*n+1:(k+1)*n);        
        u_bar_k=u_bar((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_bar_new-system_f(theta_star(1),theta_star(2),u_bar_k,x_bar_k(1),x_bar_k(2));
        con = [con; ceqnew];
    end
    %Constraints
    for k=1:N
        u_bar_k=u_bar((k-1)*m+1:k*m);
        x_bar_k=x_bar((k-1)*n+1:k*n);
        cineqnew=L*[u_bar_k;x_bar_k]-l;
        con = [con; cineqnew];
    end
    %Terminal Constraints
    x_bar_N=x_bar(end-n+1:end);
    ceqnew=x_bar_N-x_s;
    con = [con; ceqnew];
end


function cost = costfunction(N, y, x_s, u_s, Q, R,n,m)
    % Formulate the cost function to be minimized
    cost = 0;
    x_bar=y(1:n*(N+1));
    u_bar=y(n*(N+1)+1:n*(N+1)+m*N);
    
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        %controller u_hat=K*(x_bar-x_hat)+u_bar
        u_bar_k=u_bar((k-1)*m+1:k*m);
        x_bar_k=x_bar((k-1)*n+1:k*n);
        %running kost at k
        cost = cost + runningcosts(x_bar_k, u_bar_k, x_s, u_s, Q, R);
    end
    
end

function cost = runningcosts(x, u, x_eq, u_eq, Q, R)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
end


function x_kplus=simulation(x_k,u_k,theta_star)
    x_kplus=system_f(theta_star(1),theta_star(2),u_k,x_k(1),x_k(2));
end