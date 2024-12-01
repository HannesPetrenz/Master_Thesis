function [X,U,S,J,X_hat_OL,X_bar_OL,V_OL,S_OL,J_wc,time,Theta_HC_t]=RALMPC_Online(x0,SS,Q_func,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,K,F,G,d_bar,L_B,c_max,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,Ts,numberitertions)
options = optimset('Display','none',...
    'TolFun', 1e-8,...
    'MaxIter', 10000,...
    'TolConSQP', 1e-6);
%Compute the neceassary system matricis
A_cl_0=A_0+B_0*K;
A_cl_theta1=A_1+B_1*K;
A_cl_theta2=A_2+B_2*K;
%Build the constraints for y=[x_bar;x_hat;v;s;s_tilde]
%size=n*(N+1)+n*(N+1)+m*N+(N+1)+1
[A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,c_max,L_B,d_bar,H,B_p,n,m,N,p);
%Cost H
H_cost=getcost(Q,R,K,n,N,m);
%initial guess
y_init=[zeros(2*(N+1)*n+N*m,1)];
for k=0:N-1
   y_init=[y_init;(1-(rho_theta0+eta_0*L_B)^k)/(1-(rho_theta0+eta_0*L_B))*d_bar];
end

%Define the requiered variables
Theta_HC_t=Theta_HC0;
eta_t=eta_0;
theta_bar_t=theta_bar0;
theta_hat_t=theta_bar0;
rho_theta_t=rho_theta0;

for h=1:numberitertions
    tic
    t=0;
    time_array=[];
    X_array=[];
    U_array=[];
    S_array=[];
    J_array=0;
    xmeasure = x0;
    terminate=true;
    Q_wc=inf;
    while t<20
        %Set membership estimation and Point estimation
        if t>0
           %Update Theta_HC_t
           [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,H_w,h_w,A_0,A_1,A_2,B_0,B_1,B_2,p,B_p,options); 
           %Update theta_hat_t
           theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,A_0,A_1,A_2,B_0,B_1,B_2,mu,p,options);
           %Update rho_theta_t
           rho_theta_t=get_updaterho(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2,K,H,options);
        end
        [x_bar_OL,x_hat_OL,v_OL,s_OL,Q_wc]=solve_RALMPC(xmeasure,SS,Q_func,H_cost,y_init,n,m,N,theta_bar_t,eta_t,rho_theta_t,theta_hat_t,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF,Q_wc,H,options);
        %Store the data
        time_array=[time_array;t];
        X_array=[X_array,xmeasure];
        U_array=[U_array,v_OL(1,:)+K*xmeasure];
        S_array=[S_array,s_OL(1,:)];
        %open loop for sample set
        X_hat_OL_cell{t+1}=x_hat_OL;
        X_bar_OL_cell{t+1}=x_bar_OL;
        V_OL_cell{t+1}=v_OL;
        S_OL_cell{t+1}=s_OL;
        J_wc_cell{t+1}=get_worstcasecosttogo(x_hat_OL,s_OL,v_OL,H,Q_wc,Q,R,K);
        x_tminus=xmeasure;
        u_tminus=v_OL(:,1)+K*xmeasure;
        %Simulate the uncertain system
        xmeasure=dynamic(x_tminus,u_tminus,A_star,B_star,W_V);
        %Compute cost function
        J_array=J_array+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
        %Check termination condition
        terminate=check_terminalcondition(X_array(:,end),S_array(end),c_max,H);
        %Update time
        disp(t)
        t=t+1;
        checkcost=0;
        for i=1:7
            u=K*x_bar_OL(:,i)+v_OL(:,i);
            checkcost=checkcost+x_bar_OL(:,i)'*Q*x_bar_OL(:,i)+u*R*u;
        end
    end 
    %Compute the real time
    time_array=(time_array-1)./Ts;
    %close loop
    time{h}=time_array;
    X{h}=X_array;
    U{h}=U_array;
    S{h}=S_array;
    J{h}=J_array;
    %open loop for sample set
    X_hat_OL{h}=X_hat_OL_cell;
    X_bar_OL{h}=X_bar_OL_cell;
    V_OL{h}=V_OL_cell;
    S_OL{h}=S_OL_cell;
    J_wc{h}=J_wc_cell;
    %Update the sample set and the function Q
    [SS,Q_func]=get_updateSampleSet(SS,Q_func,X_bar_OL{h},J_wc{h},V_OL{h},S_OL{h});
    toc
end
end
%% Help functions
function D=get_D(x,u,A_1,A_2,B_1,B_2)
    D = [A_1*x + B_1*u, A_2*x + B_2*u];
end

function [A,B]=get_systemmatrix(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function terminate=check_terminalcondition(x,s,c_max,H)
    condition=[];
    for i=1:size(H,1)
        condition=[condition;c_max*(s+H(i,:)*x)-1];
    end
    if max(condition)<=0
        terminate=0;
    else
        terminate=1;
    end    
end

function x_tplus=dynamic(x_t,u_t,A_star,B_star,W_V)
%This function simulates the system dynamics with the disturbance
    w = W_V(1,:)';
    x_tplus=A_star*x_t+B_star*u_t+(1-2*rand)*w;
end

function [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,H_w,h_w,A_0,A_1,A_2,B_0,B_1,B_2,p,B_p,options)
    %This function updates the hypercube with Algorithm 1
    Theta_HC_tminus=Theta_HC_t;
    eta_tminus = eta_t;
    theta_bar_tminus = theta_bar_t;
    %Update Delta
    for i=1:length(Delta)-1
        Delta{i}=Delta{i+1};
    end
    %Compute the latest Delta
    h_delta = h_w - H_w*(xmeasure - A_0*x_tminus - B_0*u_tminus);
    H_delta = -H_w*get_D(x_tminus,u_tminus,A_1,A_2,B_1,B_2);
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
        eta_t = round(max(theta_max - theta_min),5);
        %Compute the new theta_bar_t
        theta_bar_t = 0.5*(theta_max + theta_min)';
        for i=1:p   %Projection    
           if theta_bar_t(i)< theta_bar_tminus(i)-0.5*(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)-0.5*(eta_tminus-eta_t);
           elseif theta_bar_t(i)>theta_bar_tminus(i)+0.5*(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)+0.5*(eta_tminus-eta_t);             
           else
             %no projection necessary  
           end          
        end
    Theta_HC_t=theta_bar_t+eta_t*B_p;
end

function theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,A_0,A_1,A_2,B_0,B_1,B_2,mu,p,options)
%This function performs the point estimate
    theta_hat_tminus=theta_hat_t;
    [A_theta_hat_t,B_theta_hat_t]=get_systemmatrix(theta_hat_tminus,A_0,A_1,A_2,B_0,B_1,B_2);
    %predicted state
    x_hat_1t=A_theta_hat_t*x_tminus+B_theta_hat_t*u_tminus;
    %Update theta_tilde
    theta_tilde_t = theta_hat_tminus + mu*get_D(x_tminus,u_tminus,A_1,A_2,B_1,B_2)'*(xmeasure-x_hat_1t);
    [theta_hat_t,~] = fmincon(@(theta) norm(theta-theta_tilde_t),zeros(p,1),Theta_HC_t.A,Theta_HC_t.b,[],[],[],[],[],options);
end

function rho_theta_t=get_updaterho(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2,K,H,options)
%This function updates the contraction rate for the updated theta_bar_t
    rho_array = [];
    [A_theta_bar_t,B_theta_bar_t] = get_systemmatrix(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2);
    A_cl_theta0 = A_theta_bar_t + B_theta_bar_t*K;
    for i = 1:size(H,1)
        [~,rho_theta_i] = linprog(-H(i,:)*A_cl_theta0, H, ones(size(H,1),1),[],[],[],[],options);
        rho_array = [rho_array; -rho_theta_i];
    end
    rho_theta_t = max(rho_array);
end

function J_wc_t=get_worstcasecosttogo(X_hat_OL,S_hat_OL,V_OL,H,Q_wc,Q,R,K)
%Compute the worst-case cost-to-go
%terminal cost
    cost_to_go_wc=Q_wc;%Set the terminal worst case cost to go. 
    %stage cost
    for k=length(X_hat_OL)-1:-1:2
        X_k_t=Polyhedron(H,S_hat_OL(k)*ones(size(H,1),1)+H*X_hat_OL(:,k));
        cost=0;
        for i=1:size(X_k_t.V,1)
            cost=max([cost;X_k_t.V(i,:)*Q*X_k_t.V(i,:)'+(X_k_t.V(i,:))*K'*R*K*X_k_t.V(i,:)'+V_OL(:,k)*R*V_OL(:,k)]);
        end
        cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
    end
    J_wc_t=flip(cost_to_go_wc);
end


function [A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,c_max,L_B,d_bar,H,B_p,n,m,N,p)
%This function builds the equality and the inequality constraints.
%Equality constraints
%System dynamics for x_bar
    A_eq=[eye(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1),zeros(n,1)]; %inital condition
    A_eq_theta1_bar=[zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    A_eq_theta2_bar=[zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    b_eq=[x0];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*(k)),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta1_bar=[A_eq_theta1_bar;zeros(n,n*k),A_cl_theta1,zeros(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*k),B_1,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta2_bar=[A_eq_theta2_bar;zeros(n,n*k),A_cl_theta2,zeros(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*k),B_2,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];    
        b_eq=[b_eq;zeros(n,1)];
    end
    %System dynamics for x_hat
    A_eq=[A_eq;zeros(n,n*(N+1)),eye(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    b_eq=[b_eq;x0];
    A_eq_theta1_hat=[zeros(size(A_eq_theta1_bar));zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    A_eq_theta2_hat=[zeros(size(A_eq_theta2_bar));zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta1_hat=[A_eq_theta1_hat;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_theta1,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_1,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta2_hat=[A_eq_theta2_hat;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_theta2,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_2,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        b_eq=[b_eq;zeros(n,1)];
    end
    %Inital condition s
    A_eq=[A_eq;zeros(1,2*n*(N+1)+N*m),1,zeros(1,N),zeros(1,1)];
    b_eq=[b_eq;0];
    A_eq_theta1_hat=[A_eq_theta1_hat;zeros(1,n*(N+1)+n*(N+1)+m*N+(N+1)),zeros(1,1)];
    A_eq_theta2_hat=[A_eq_theta2_hat;zeros(1,n*(N+1)+n*(N+1)+m*N+(N+1)),zeros(1,1)];
    A_eq_theta1_bar=[A_eq_theta1_bar;zeros(-size(A_eq_theta1_bar,1)+size(A_eq,1),size(A_eq_theta1_bar,2))];
    A_eq_theta2_bar=[A_eq_theta2_bar;zeros(-size(A_eq_theta2_bar,1)+size(A_eq,1),size(A_eq_theta2_bar,2))];
    %Inequality constraints
    %Terminal constraints
    A_ineq=[zeros(size(H,1),n*N),H,zeros(size(H,1),n*(N+1)),zeros(size(H,1),N*m),zeros(size(H,1),N+1),-ones(size(H,1),1)];
    b_ineq_xF=[H];
    A_ineq=[A_ineq;zeros(1,n*(N+1)),zeros(1,n*(N+1)),zeros(1,N*m),zeros(1,N),1,1];
    b_ineq_sF=[zeros(size(H,1),1);1];
    %Constraints
    b_ineq=zeros(size(A_ineq,1),1);
    for k=0:N-1
        A_ineq=[A_ineq;zeros(size(F,1),n*k),F+G*K,zeros(size(F,1),n*(N-k)),zeros(size(F,1),n*(N+1)),zeros(size(F,1),m*k),G,zeros(size(F,1),m*(N-1-k)),zeros(size(F,1),k),c,zeros(size(F,1),N-k),zeros(size(F,1),1)];
        b_ineq=[b_ineq;ones(size(F,1),1)];
    end
    %Tube dynamic s
    A_ineq_eta=zeros(size(A_ineq));
    A_ineq_rho=zeros(size(A_ineq));
    for k=0:N-1
        for j=1:2^p  
                A_ineq_eta=[A_ineq_eta;...
              zeros(size(H,1),n*k),H*(B_p.V(j,1)*A_cl_theta1+B_p.V(j,2)*A_cl_theta2),zeros(size(H,1),n*(N-k)),zeros(size(H,1),n*(N+1)),zeros(size(H,1),m*k),H*(B_p.V(j,1)*B_1+B_p.V(j,2)*B_2),zeros(size(H,1),m*(N-1-k)),zeros(size(H,1),k),L_B*ones(size(H,1),1),zeros(size(H,1),N-k),zeros(size(H,1),1)];
        end
        b_ineq=[b_ineq;-d_bar*ones(2^p*size(H,1),1)];
        A_ineq_rho=[A_ineq_rho;...
             zeros(2^p*size(H,1),2*n*(N+1)+m*N),zeros(2^p*size(H,1),k),ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k),zeros(2^p*size(H,1),1)];
        A_ineq=[A_ineq;...
             zeros(2^p*size(H,1),2*n*(N+1)+m*N),zeros(2^p*size(H,1),k+1),-ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k-1),zeros(2^p*size(H,1),1)];
    end
    b_ineq_xF=[b_ineq_xF;zeros(size(A_ineq,1)-size(b_ineq_xF,1),size(b_ineq_xF,2))];
    b_ineq_sF=[b_ineq_sF;zeros(size(A_ineq,1)-size(b_ineq_sF,1),size(b_ineq_sF,2))];
end

function H_cost=getcost(Q,R,K,n,N,m)
%This function builds the cost function matrix H
    H_cost=[];
    %Zeros for x_bar
    H_cost=blkdiag(H_cost,zeros(n*(N+1)));
    %Cost for x_hat
    for k=0:N-1
       H_cost=blkdiag(H_cost,Q+K'*R*K); 
    end
    %Terminal cost
    H_cost=blkdiag(H_cost,zeros(size(Q)));
    %Cost for v
    for k=0:N-1
       H_cost=blkdiag(H_cost,R); 
    end
    %Add cross-terms
    for k=0:N-1
        H_cost(2*(N+1)*n+k+m,(N+1)*n+k*n+1:(N+1)*n+(k+1)*n) = R*K;
        H_cost((N+1)*n+k*n+1:(N+1)*n+(k+1)*n,2*(N+1)*n+k+1) = K'*R;
    end
    %Zeros for s
    H_cost=blkdiag(H_cost,zeros(N+1+1));
    H_cost=2*H_cost; %Since we solve the problem min 0.5*y'*H*y
end

function [SS,Q_func]=get_updateSampleSet(SS,Q_func,X_bar_cell,J_wc_cell,V_OL_cell,S_OL_cell)
    %This function updates the Sample set and Q_fun
    for i=1:length(X_bar_cell)
        SS=[SS,[X_bar_cell{1}(:,2:end-1);S_OL_cell{1}(2:end-1)';V_OL_cell{1}(:,2:end)]];
        Q_func=[Q_func,J_wc_cell{1}(2:end)];
    end
end