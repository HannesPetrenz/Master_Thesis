function [X,U,J,X_hat_OL,X_bar_OL,U_bar_OL,S_OL,J_wc,Theta_HC,theta_bar,theta_hat,eta,t_cpu]=RALMPC(x0,x_s, u_s, Q, R,P,Theta_HC0,theta_star,theta_bar0,eta_0,rho_theta0,L_B_rho,d_bar,L,l,Delta,D,options,mu,B_p,p,r,N,numberitertions,mpciterations,m,n,c,SS,Q_func)
import casadi.*
%Option for casadi solver
opts = struct;
opts.ipopt.print_level = 0;  % Suppress solver output
opts.print_time = false;      % Disable timing information
% make symbolic x_bar x_hat u_bar s w s_tilde (n*(N+1)+n*(N+1)+m*N+(N+1)+2^p*N+1)
y=MX.sym('y',n*(N+1)+n*(N+1)+m*N+(N+1)+N+1);
%nonlinear constraints 
con_lb=[zeros(n*(N),1);zeros(n*(N),1);zeros((N),1);zeros((N),1);-inf*ones(r*N,1);-inf*ones(2,1)];
con_ub=[zeros(n*(N),1);zeros(n*(N),1);zeros((N),1);zeros((N),1);zeros(r*N,1);zeros(2,1)];
%input constraints: use lb, ub
lb=[x0;-inf*ones(n*(N),1);x0;-inf*ones(n*(N),1);-inf*ones(m*(N),1);zeros(N+1,1);zeros(N,1);zeros(1,1)];
ub=[x0;inf*ones(n*(N),1);x0;inf*ones(n*(N),1);inf*ones(m*(N),1);0;inf*ones((N),1);inf*ones((N),1);inf*ones(1,1)];
%Define the requiered variables
Theta_HC_t=Theta_HC0;
eta_t=eta_0;
theta_bar_t=theta_bar0;
theta_hat_t=theta_bar0;
rho_theta_t=rho_theta0;

%Iteration over Learning iteration
for h=1:numberitertions
    t=0;
    J{h}{1}=0;
    xmeasure = x0;
    tStart = cputime;
    X{h}=[];
    U{h}=[];
    %Initial solution
    L_Theta=eta_t*L_B_rho;
    y_init=[repmat(x0,N+1,1);repmat(x0,N+1,1);zeros(m*N,1)];
    for k=0:N-1
       y_init=[y_init;(1-(rho_theta0+L_Theta)^k)/(1-(rho_theta0+L_Theta))*d_bar+10^-4];
    end
    y_init=[y_init;0;d_bar*ones(N,1);0];
    
        while t<mpciterations %MPC iteration
            idx_infeasible=[];
            
            %Set membership estimation and Point estimation
            if t>0 
               [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,D,options,B_p,p);
               theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,mu,p,options);
               rho_theta_t=rho_theta0+(eta_0-eta_t)*L_B_rho;
            end
            % Set initial constraint
            lb(1:n)=xmeasure;
            ub(1:n)=xmeasure;
            lb(n*(N+1)+1:n*(N+1)+n)=xmeasure;
            ub(n*(N+1)+1:n*(N+1)+n)=xmeasure;
            %pre-select only the promising values in SS and Q_func
            if t>0
                if h>1
                    SS_select=[SS{1},SS{end}]; %only the iteration 0 and h-1
                    Q_func_select=[Q_func{1},Q_func{end}];%only the iteration 0 and h-1
                else
                    SS_select=[SS{1}]; %only the iteration 0 and h-1
                    Q_func_select=[Q_func{1}];%only the iteration 0 and h-1
                end
                idx=find(Q_func_select<=Q_current); %cost to go must be larger equl to the previous selected
                SS_select=SS_select(:,idx);
                Q_func_select=Q_func_select(:,idx);
                y_best=[];
                cost_best=inf;
                Q_current=inf;
                SS_current=[];
            else
                if h>1
                    SS_select=[SS{1},SS{end}]; %only the iteration 0 and h-1
                    Q_func_select=[Q_func{1},Q_func{end}];%only the iteration 0 and h-1
                else
                    SS_select=[SS{1}]; %only the iteration 0 and h-1
                    Q_func_select=[Q_func{1}];%only the iteration 0 and h-1
                end
                y_best=[];
                cost_best=inf;
                Q_current=inf;
                SS_current=[];
            end
            
            %Iteration over all values in SS_selct to determine the best
            %terminal value
            for ii=1:size(SS_select,2)
                if ~ismember(ii, idx_infeasible) %Check if the candidate in SS_select is in the infeasible region
                    % Update the constrain set and the optimization problem
                    [con] = nonlinearconstraints(N, y, n, m,theta_bar_t,theta_hat_t,eta_t,L_B_rho,rho_theta_t,L,l,c,SS_select(1:n,ii),SS_select(n+1,ii),P,d_bar);
                    % objective function
                    obj=costfunction(N, y, x_s, u_s, Q, R,n,Q_func_select(ii),m);
                    nlp = struct('x', y, 'f', obj, 'g', con);
                    %construct optimization problem
                    solver = nlpsol('solver', 'ipopt', nlp, opts);
                    %Solve the optimization problem
                    res = solver('x0' , y_init,... % solution guess
                             'lbx', lb,...           % lower bound on x
                             'ubx', ub,...           % upper bound on x
                             'lbg', con_lb,...           % lower bound on g
                             'ubg', con_ub);             % upper bound on g
                    %check if the solver was succesful
                    stats = solver.stats();
                    %Check if the optimization problem is succesfully
                    %solved 
                    if strcmp(stats.return_status, 'Solve_Succeeded')
                        %disp('Optimization was successful!');
                        %Check if the current vlaue of SS is better then
                        %the best previous one and update in case
                        y_current=full(res.x);
                        cost_current=costfunction(N, y_current, x_s, u_s, Q, R,n,Q_func_select(ii),m);
                        if cost_best>= cost_current
                            y_best=y_current;
                            cost_best=cost_current;
                            Q_current=Q_func_select(ii);
                            SS_current=SS{h}(:,ii);
                        end
                    elseif strcmp(stats.return_status, 'Infeasible_Problem_Detected')
                            %disp(['Solver failed with status: ', stats.return_status]);
                            %Check if the value is infeasible and update
                            %infeasible list
                            % if ~ismember(ii, idx_infeasible) && ~isempty(idx_infeasible)
                            %     disp("shit")
                            % end
                            if isempty(idx_infeasible)
                                %idx_infeasible=find(norm(SS_select(1:n,ii)-x_s)>=vecnorm(SS_select(1:n,:)-x_s));
                            end
                    elseif strcmp(stats.return_status, 'Maximum_Iterations_Exceeded')
                           %disp(['Solver failed with status: ', stats.return_status]);
                    else
                           disp(['Solver failed with status: ', stats.return_status]);
                    end
                end
            end


            %Open loop solution
            y_OL=y_best; 
            X_bar_OL{h}{t+1}=reshape(y_OL(1:n*(N+1)),2,[]);
            X_hat_OL{h}{t+1}=reshape(y_OL(n*(N+1)+1:n*(N+1)+n*(N+1)),2,[]);
            U_bar_OL{h}{t+1}=y_OL(n*(N+1)+n*(N+1)+1:n*(N+1)+n*(N+1)+m*N);
            S_OL{h}{t+1}=y_OL(n*(N+1)+n*(N+1)+m*N+1:n*(N+1)+n*(N+1)+m*N+N+1);
            J_wc{h}{t+1}=get_worstcasecosttogo(X_hat_OL{h}{1},S_OL{h}{1},U_bar_OL{h}{1},P,Q_current,Q,R,L(:,m+1:end),l);
            %plotting(X_bar_OL{h}{t+1},S_OL{h}{t+1},P,SS,n,SS_current,y_OL(end),Q_current)
            %Parameter update 
            Theta_HC{h}{t+1}=Theta_HC_t;
            theta_bar{h}{t+1}=theta_bar_t;
            theta_hat{h}{t+1}=theta_hat_t;
            eta{h}{t+1}=eta_t;
            %Close Loop
            u_cl=U_bar_OL{h}{end}(1);
            X{h}=[X{h},xmeasure];
            U{h}=[U{h},u_cl];
            x_tminus=xmeasure;
            u_tminus=u_cl;
            xmeasure=simulation(X{h}(:,end),u_cl,theta_star,D,n);
            %New initial solution for the optimization
            u_init = controller_K(u_s,x_s(1),x_s(2))*(X_bar_OL{h}{end}(:,end)-x_s)+u_s;
            x_init = system_f(theta_bar_t(1),theta_bar_t(2),u_init,X_bar_OL{h}{end}(1,end),X_bar_OL{h}{end}(2,end));
            s_init=rho_theta_t*S_OL{h}{end}(end)+eta_t*L_B_rho*S_OL{h}{end}(end)+d_bar+max(uncertainty_w_Theta(eta_t,X_bar_OL{h}{end}(1,end),X_bar_OL{h}{end}(2,end)));
            w_init=max(uncertainty_w_Theta(eta_t,X_bar_OL{h}{end}(1,end),X_bar_OL{h}{end}(2,end)));
            y_init=[y_OL(n+1:n*(N+1));x_init;y_OL(n*(N+1)+n+1:n*(N+1)+n*(N+1));x_init;y_OL(n*(N+1)+n*(N+1)+m+1:n*(N+1)+n*(N+1)+m*N);u_init;y_OL(n*(N+1)+n*(N+1)+m*N+1+1:n*(N+1)+n*(N+1)+m*N+N+1);s_init;y_OL(n*(N+1)+n*(N+1)+m*N+N+1+1+1:end);w_init];
            %Compute cost function
            J{h}{end+1}=J{h}{end}+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
            %Next step
            t=t+1;
        end
        %Update the sample set and the function Q
    [SS,Q_func]=get_updateSampleSet(SS,Q_func,X_bar_OL{h},J_wc{h},S_OL{h},U_bar_OL{h},h);
    t_cpu{h} = cputime - tStart;
end
end


function [con] = nonlinearconstraints(N, y, n, m,theta_bar,theta_hat,eta_t,L_B_rho,rho_thetat,L,l,c,x_s,s_s,P,d_bar)
    % Introduce the equality constraints also for the terminal state  
       x_bar=y(1:n*(N+1));
       x_hat=y(n*(N+1)+1:n*(N+1)+n*(N+1));
       u_bar=y(n*(N+1)+n*(N+1)+1:n*(N+1)+n*(N+1)+m*N);
       s=y(n*(N+1)+n*(N+1)+m*N+1:n*(N+1)+n*(N+1)+m*N+N+1);
       w=y(n*(N+1)+n*(N+1)+m*N+N+1+1:end-1);
       s_tilde=y(end);
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
    
    cineqnew=(x_bar_N-x_s).'*P*(x_bar_N-x_s)-s_tilde^2;
    con = [con; cineqnew];
    cineqnew=s_N-s_s+s_tilde;
    con = [con; cineqnew];
end


function cost = costfunction(N, y, x_s, u_s, Q, R,n,Q_func,m)
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
    cost = cost + Q_func;
    
end

function cost = runningcosts(x, u, x_eq, u_eq, Q, R)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
end


function [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,D,options,B_p,p)
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

function theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,mu,p,options)
%This function performs the point estimate
    theta_hat_tminus=theta_hat_t;
    %predicted state
    x_hat_1t=system_f(theta_hat_tminus(1),theta_hat_tminus(2),u_tminus,x_tminus(1),x_tminus(2));
    %Update theta_tilde
    theta_tilde_t = theta_hat_tminus + mu*system_G(x_tminus(1),x_tminus(2))*(xmeasure-x_hat_1t);
    [theta_hat_t,~] = fmincon(@(theta) norm(theta-theta_tilde_t),zeros(p,1),Theta_HC_t.A,Theta_HC_t.b,[],[],[],[],[],options);
end

function [SS,Q_func]=get_updateSampleSet(SS,Q_func,X_bar_OL_h,J_wc_h,S_OL_h,U_bar_OL,h)
    %This function updates the Sample set and Q_fun
    SS{h+1}=[];
    Q_func{h+1}=[];
    for i=1:length(X_bar_OL_h)
        SS{h+1}=[SS{h+1},[X_bar_OL_h{i}(:,2:end-1);S_OL_h{i}(2:end-1,:)';U_bar_OL{i}(2:end,:)']];
        Q_func{h+1}=[Q_func{h+1},J_wc_h{i}];
    end
end

function J_wc_t=get_worstcasecosttogo(X_hat_OL,S_hat_OL,U_bar_OL,P,Q_wc,Q,R,L_x,l_x)
%Compute the worst-case cost-to-go
%terminal cost
    cost_to_go_wc=Q_wc;%Set the terminal worst case cost to go. 
    %stage cost
    for k=length(X_hat_OL)-1:-1:2
        K = controller_K(U_bar_OL(k,:),X_hat_OL(1,k),X_hat_OL(2,k));
        cost=worstcost(X_hat_OL(:,k),U_bar_OL(k,:),P,S_hat_OL(k,:),Q,R,K,L_x,l_x);
        cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
    end
    J_wc_t=flip(cost_to_go_wc(2:end));
end

function J_wc=worstcost(z,v,P,s,Q,R,K,L_x,l_x)

options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
x = fmincon(@(x)costfunction(x,v,Q,R,K),z,L_x,l_x,[],[],[],[],@(x)nonlinearcon(x,P,z,s), options);
J_wc=-(costfunction(x,v,Q,R,K))+v'*R*v;

    function cost=costfunction(x,v,Q,R,K)
        Q_bar=Q+K'*R*K;
        c=2*v'*R*K;
        cost=-(x'*Q_bar*x+c*x);
    end
    function [c,ceq] =nonlinearcon(x,P,z,s)
        if cond(P) > 1e6
            P = P + 1e-6 * eye(size(P)); % Regularization for numerical stability
        end
    c=(x-z)'*P*(x-z)-s^2;
    ceq = [];
    end
end

function x_kplus=simulation(x_k,u_k,theta_star,D,n)
    w = D.V(1,:)';
    disturbance=(1-2*rand(n,1)).*w;
    x_kplus=system_f(theta_star(1),theta_star(2),u_k,x_k(1),x_k(2))+disturbance;
end


function plotting(X_OL,S_OL,P,SS,n,SS_current,s_tilde,Q)
figure;
hold on
Polyhedron_Color = [0.8, 0.8, 1];  % Light blue for polyhedron (excluded from legend)
Polyhedron_Color_2 = [1, 0.6, 0.2];  % Light blue for polyhedron (excluded from legend)
X_Init_Color = [0.2, 0.2, 0.8];    % Dark blue for initial trajectory
RALMPC_Color = [0.6, 0.8, 0.2]; % Light Green (RALMPC First Iteration)
SS=SS{1};

% Plot the Polyhedron constraints (EXCLUDED from the legend)
numberpoints=100;
for k = 1:length(SS)
    X_ellipse = pointsellipse(P, SS(1:n, k), SS(n+1,k), numberpoints);
    h = patch(X_ellipse(1, :), X_ellipse(2, :), Polyhedron_Color, ...
        'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 2, ...
        'FaceAlpha', 0.3); % Adjust transparency
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); % Exclude from legend
end

% Plot initial trajectory (X_bar_initial)
plot(SS(1,:), SS(2,:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'Color', X_Init_Color, 'DisplayName', 'Initial Trajectory');

% Plot the Polyhedron constraints (EXCLUDED from the legend)
numberpoints=100;
for k = 1:length(X_OL)
    X_ellipse = pointsellipse(P, X_OL(:, k), S_OL(k), numberpoints);
    h = patch(X_ellipse(1, :), X_ellipse(2, :), Polyhedron_Color_2, ...
        'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5, ...
        'FaceAlpha', 0.3); % Adjust transparency
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); % Exclude from legend
end
% Plot RAMPC trajectory
plot(X_OL(1, :), X_OL(2, :), 'Marker', '.','LineStyle','--', 'LineWidth', 2, 'DisplayName', 'RAMPC', 'Color', RALMPC_Color);
grid on
(X_OL(:,end)-SS_current(1:2))'*P*(X_OL(:,end)-SS_current(1:2))-s_tilde
S_OL(end)-SS_current(3)+s_tilde
function ellipse=pointsellipse(P,center,delta,numberpoints)
[eig_U, eig_D] = eig(P);
theta_circle = linspace(0, 2*pi, numberpoints);
circle = [cos(theta_circle); sin(theta_circle)];
A_trans = eig_U * sqrt(inv(eig_D))*delta;
ellipse = A_trans * circle + center;
end
end