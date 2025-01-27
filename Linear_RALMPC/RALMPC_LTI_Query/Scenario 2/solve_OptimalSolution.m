function [X_bar_OL,V_OL,J,X,U,time]=solve_OptimalSolution(x0,A_star,B_star,K,F,G,Q,R,n,m,N,W_V)
options = optimset('Display','none',...
    'TolFun', 1e-8,...
    'MaxIter', 10000,...
    'TolConSQP', 1e-6);
%Compute the neceassary system matricis
A_cl=A_star+B_star*K;
%Build the constraints for y=[x_bar;v;s] size=n*(N+1)+m*N
[A_eq,b_eq,A_ineq,b_ineq]=get_constraints(x0,A_cl,B_star,K,F,G,n,m,N);
%Cost H
H_cost=getcost(Q,R,K,n,N,m);
%MPC Iteration
%initial guess
y_init=[zeros((N+1)*n+N*m,1)];

%MPC Iteration
%initial guess
t=0;
xmeasure = x0;
time=[];
terminate=true;
X=[];
U=[];
S=[];
J=0;
%MPC iteration loop
while t<60
    %Optimization
    %Update the constraints with the updated parameter
    b_eq(1:n)=xmeasure;
    %Optimize the quadratic programm 
    [y_OL,V,exitflag] = quadprog(H_cost,zeros(n*(N+1)+m*N,1),A_ineq,b_ineq,A_eq,b_eq,[],[],y_init,options);
    %Extract the open loop solution
    x_bar_OL = reshape(y_OL(1:2*(N+1))',n,[]);
    v_OL=reshape(y_OL((N+1)*n+1:(N+1)*n+m*N),m,[]);
    %Store the data
        %close loop
        time=[time;t];
        X=[X,xmeasure];
        U=[U,v_OL(1:m)+K*xmeasure];
        %open loop for sample set
        X_bar_OL{t+1}=x_bar_OL;
        V_OL{t+1}=v_OL;
    x_tminus=X(:,end);
    u_tminus=U(:,end);
    J=J+x_tminus'*Q*x_tminus+x_tminus'*R*x_tminus;
    %Simulate the uncertain system
    xmeasure=dynamic(X(:,end),U(end),A_star,B_star,W_V);
    t=t+1;
end
end
%% Help functions

function x_tplus=dynamic(x_t,u_t,A_star,B_star,W_V)
%This function simulates the system dynamics with the disturbance
    w = 0;%W_V(1,:)';
    x_tplus=A_star*x_t+B_star*u_t+(1-2*rand)*w;
end

function [A_eq,b_eq,A_ineq,b_ineq]=get_constraints(x0,A_cl,B_star,K,F,G,n,m,N)
%This function builds the equality and the inequality constraints.
%Equality constraints
%System dynamics for x_bar
    A_eq=[eye(n),zeros(n,n*N),zeros(n,m*N)]; %inital condition
    b_eq=[x0];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*k),A_cl,-eye(n),zeros(n,n*(N-k-1)),zeros(n,m*(k)),B_star,zeros(n,m*(N-k-1))];
        b_eq=[b_eq;zeros(n,1)];
    end
    
    %Inequality constraints
    %Terminal constraints
    A_ineq=[];
    b_ineq=[];
    %Constraints
    for k=0:N-1
        A_ineq=[A_ineq;zeros(size(F,1),n*k),F+G*K,zeros(size(F,1),n*(N-k)),zeros(size(F,1),m*k),G,zeros(size(F,1),m*(N-1-k))];
        b_ineq=[b_ineq;ones(size(F,1),1)];
    end
end

function H_cost=getcost(Q,R,K,n,N,m)
%This function builds the cost function matrix H
    H_cost=[];
    %Cost for x_hat
    for k=0:N
       H_cost=blkdiag(H_cost,Q+K'*R*K); 
    end
   
    for k=0:N-1
       H_cost=blkdiag(H_cost,R); 
    end
    %Add cross-terms
    for k=0:N-1
        H_cost(1*(N+1)*n+k+m,k*n+1:+(k+1)*n) = R*K;
        H_cost(k*n+1:(k+1)*n,(N+1)*n+k+1) = K'*R;
    end
    H_cost=2*H_cost; %Since we solve the problem min 0.5*y'*H*y
end