function ARMPC
% Linear system, 3 constant parameters, RMPC with parameter estimation, ML
% example. X_0 PI lambda-contractive, X_T RPI, RLS estimation, simpler
% alpha evolution 
% ONLINE COMPUTATIONS
clear all
close all
clc
%% Load and import
load('ARMPC_Offline.mat');
disp('Offline things loaded');
disp('Quadprog')
%error('make correct cost with u=v+k*x')
adapt=true;
disp(['Adaptation: ' num2str(adapt)])
tol_opt       = 1e-8;
options = optimset('Display','none',...
    'TolFun', tol_opt,...
    'MaxIter', 10000,...
    'TolConSQP', 1e-6);
%    'Algorithm', 'active-set',...
t=tic;
% just some variables I need
options_optimize = sdpsettings('verbose',0);
T_lim = abs(Theta.V(1,1));  % I need this to plot Theta^HC_t
not_contained_in_prev_HC = 0;   % to check if Theta^HC_t cintained in Theta^HC_t-1

Theta_HC = Theta;
HC_A = Theta_HC.A;
HC_b = Theta_HC.b;
%% MPC things
N = 14;
M = 10;
step_width = 25;                    % set reference step width
MPCiterations = step_width*5;
x_init = [0;0];                     % Initial state
t_init = 0.0;                       % Initial time

% steady state
x_ss = [1;0];
u_ss = x_ss(1);

%% Parameter init
% side and center of Theta
%hypercube
eta=2; %start with \|\theta\|_\infty\leq 1
theta_overline = zeros(p,1);
%LMS
%theta_hat = theta_overline;
theta_RLS = theta_overline;
HC_center = theta_overline;
HC_side = eta;
rho_t=compute_rho(theta_overline,A_0,A_1,A_2,B_0,B_1,B_2,Hx,K);
% one of the verteces oh Theta to have real theta almost in the corner
theta_true = [1;-1];%zeros(p,1);%Theta.V(1,:)'*0.8;
% I fill M Delta sets with Theta (for t < M)
for i = 0:(M-1)  
    eval(['Delta_' num2str(i) '= Theta;']);
end
theta_out = 0;
%% Constraints
disp('Building linear constraints....');
%optimization varaible y
%y=[x,\hat{x},v,s]: n*(N+1)+N*m+N+1 (add x_hat later...)
lb=[];%[repmat([x_1min;x_2min],2*(N+1),1);-inf(N*m+N+1,1)];
ub=[];%repmat([x_1max;x_2max],2*(N+1),1);inf(N*m+N+1,1)];
%A_eq_t=[A_eq...
%        A_eq_theta_hat_1*theta_hat_1+A_eq_theta_hat_2*theta_hat_2+...
%        A_eq_theta_overline_1*theta_overline_1+A_eq_theta_overline_2*theta_overline_2+...
%initial state
y_length=2*n*(N+1)+N*m+N+1;
A_eq=[eye(n),zeros(n,n*N+n*(N+1)+N*m+N+1,1);...
      zeros(n,n*(N+1)),eye(n),zeros(n,N*n+N*m+N+1)];
b_eq=[x_init;x_init];%zeros(n,1);
%initial tube
A_eq=[A_eq;zeros(1,2*n*(N+1)+N*m),1,zeros(1,N)];
b_eq=[b_eq;0];
A_eq_theta_overline_1=zeros(2*n+1,y_length);
A_eq_theta_overline_2=zeros(2*n+1,y_length);
A_eq_theta_hat_1=zeros(2*n+1,y_length);
A_eq_theta_hat_2=zeros(2*n+1,y_length);
%dynamic
%center tube -  state dynamics
A_cl_0=A_0+B_0*K;
A_cl_1=A_1+B_1*K;
A_cl_2=A_2+B_2*K;
for k=0:N-1
A_eq=[A_eq;...
    zeros(n,k*n),A_cl_0,-eye(n),zeros(n,(N-k-1)*n),zeros(n,n*(N+1)),zeros(n,k*m),B_0,zeros(n,(N-k-1)*m),zeros(n,N+1)];
A_eq_theta_overline_1=[A_eq_theta_overline_1;...
     zeros(n,k*n),A_cl_1,zeros(n,(N-k)*n),zeros(n,n*(N+1)),zeros(n,k*m),B_1,zeros(n,(N-k-1)*m),zeros(n,N+1)];
A_eq_theta_overline_2=[A_eq_theta_overline_2;...
     zeros(n,k*n),A_cl_2,zeros(n,(N-k)*n),zeros(n,n*(N+1)),zeros(n,k*m),B_2,zeros(n,(N-k-1)*m),zeros(n,N+1)];
end
b_eq=[b_eq;zeros(N*n,1)];
A_eq_theta_hat_1=[A_eq_theta_hat_1;zeros(N*n,y_length)];
A_eq_theta_hat_2=[A_eq_theta_hat_2;zeros(N*n,y_length)];
%
%dynamics LMS pint estimate
for k=0:N-1 
A_eq=[A_eq;...
    zeros(n,n*(N+1)),zeros(n,k*n),A_cl_0,-eye(n),zeros(n,(N-k-1)*n),zeros(n,k*m),B_0,zeros(n,(N-k-1)*m),zeros(n,N+1)];
A_eq_theta_hat_1=[A_eq_theta_hat_1;...
    zeros(n,n*(N+1)),zeros(n,k*n),A_cl_1,zeros(n,(N-k)*n),zeros(n,k*m),B_1,zeros(n,(N-k-1)*m),zeros(n,N+1)];
A_eq_theta_hat_2=[A_eq_theta_hat_2;...
    zeros(n,n*(N+1)),zeros(n,k*n),A_cl_2,zeros(n,(N-k)*n),zeros(n,k*m),B_2,zeros(n,(N-k-1)*m),zeros(n,N+1)];
end
b_eq=[b_eq;zeros(N*n,1)];
A_eq_theta_overline_1=[A_eq_theta_overline_1;zeros(N*n,y_length)];
A_eq_theta_overline_2=[A_eq_theta_overline_2;zeros(N*n,y_length)];
%%
%terminal constraint
A_ineq=[zeros(c_0,n*N),cmax*X_0.A,zeros(c_0,n*(N+1)+N*m),zeros(c_0,N),cmax*ones(c_0,1)];
b_ineq=[ones(c_0,1)];
b_ineq_xss=cmax*X_0.A;
%inequality constraints (tightened)
for k=0:N-1
   A_ineq=[A_ineq;...
          zeros(q,n*k),F+G*K,zeros(q,n*(N-k)),zeros(q,n*(N+1)),zeros(q,m*k),G,zeros(q,m*(N-1-k)),zeros(q,k),c_j,zeros(q,N-k)]; 
end
b_ineq=[b_ineq;ones(q*N,1)];
b_ineq_xss=[b_ineq_xss;zeros(q*N,n)];
A_ineq_eta=zeros(size(A_ineq));
A_ineq_rho=zeros(size(A_ineq));
%%
%dynamic tube size s
%define vertices e_j
B_p=Polyhedron([1,0;-1,0;0,1;0,-1],0.5*ones(2*p,1));
e_j=B_p.V';
%
for k=0:N-1
    for j=1:2^p        
    A_ineq_eta=[A_ineq_eta;...
          zeros(c_0,n*k),X_0.A*(e_j(1,j)*A_cl_1+e_j(2,j)*A_cl_2),zeros(c_0,n*(N-k)),zeros(c_0,n*(N+1)),zeros(c_0,m*k),X_0.A*(e_j(1,j)*B_1+e_j(2,j)*B_2),zeros(c_0,m*(N-1-k)),zeros(c_0,k),L_B*ones(c_0,1),zeros(c_0,N-k)];
    end
    b_ineq=[b_ineq;-dbar*ones(2^p*c_0,1)];
    A_ineq_rho=[A_ineq_rho;...
         zeros(2^p*c_0,2*n*(N+1)+m*N),zeros(2^p*c_0,k),ones(2^p*c_0,1),zeros(2^p*c_0,N-k)];
    A_ineq=[A_ineq;...
         zeros(2^p*c_0,2*n*(N+1)+m*N),zeros(2^p*c_0,k+1),-ones(2^p*c_0,1),zeros(2^p*c_0,N-k-1)];
end
b_ineq_xss=[b_ineq_xss;zeros(2^p*c_0*N,n)];

%A_ineq=[A_ineq;zeros(c_0*N*2^p,y_length)];
toc(t)
%% cost function 

H=[];
%cost robust
H=blkdiag(H,zeros(n*(N+1)));
Q_star=Q+K'*R*K;
%cost LMS
for k=0:N-1
   H=blkdiag(H,Q_star); 
end
H=blkdiag(H,P);
%cost v
for k=0:N-1
   H=blkdiag(H,R); 
end
%add cross terms: v'*R*K*x
for k=0:N-1
   H(2*(N+1)*n+k,(N+1)*n+k*n+1:(N+1)*n+(k+1)*n) = R*K;
   H((N+1)*n+k*n+1:(N+1)*n+(k+1)*n,2*(N+1)*n+k) = K'*R;
end
%cost s
H=blkdiag(H,zeros(N+1));
H=0.5*H;
%quadprog: min 0.5y'Hy+f'y
%min 0.5\|y-y_s\|_H^2=0.5y'*H*y-y'H*y_s
y_s=[repmat(x_ss,2*(N+1),1);repmat(u_ss-K*x_ss,N,1);zeros(N+1,1)];
f=-y_s'*H;
%% initial guess
y_init=[zeros(2*(N+1)*n+N*m,1)];
for k=0:N-1
   y_init=[y_init;(1-(rho_t+eta*L_B)^k)/(1-(rho_t+eta*L_B))*dbar];
end
%% Set variables for output
t = [];
x = [];
s = [];
u = [];
x1 = [];
x1_ref = [];

theta_HC_error = [];
theta_RLS_error = [];
error_OL = [];
theta_RLS_rec=[];
theta_LMS_store=[ ];
theta_center_store=[ ];
eta_store=[ ];
%% Initilization of measured values and others
tmeasure = t_init;
xmeasure = x_init;
%% SIMULATION
disp('Starting MPC, total number of iterations to make:');
disp(MPCiterations);
MPCstart = tic;
for ii = 1:MPCiterations
    %% Set point
    if ii > step_width && ii <= step_width*2
        %error('program setpoint change') 
        x_ss = [0;0];
    elseif ii > step_width*2 && ii <= step_width*3
        x_ss = [1;0];
    elseif ii > step_width*3 && ii <= step_width*4
        x_ss = [0;0];
    elseif ii > step_width*4 && ii <= step_width*5
        x_ss = [1;0];
    elseif ii > step_width*5 && ii <= step_width*6
        x_ss = [0;0];
    elseif ii > step_width*6 && ii <= step_width*7
        x_ss = [1;0];
    end
    u_ss = x_ss(1);
    %% Find Theta_t  
    if adapt==true    
     % init
    Theta_HC_prev = Theta_HC;
    %H_HC_prev = Theta_HC_prev.A;
    %h_HC_prev = Theta_HC_prev.b;
    HC_side_prev = HC_side;
    HC_center_prev = HC_center;
    
    
   % Find $\Delta_t$
%    error('programm parameter update (needs to be slightly adjusted with eta..B_p)')
    if ii>1
        % I must update the delta sets (shift them all)
        for i = 0:(M-2) 
            eval(['Delta_' num2str(M-1-i) '= Delta_' num2str(M-2-i) ';']);
        end
        h_delta = h_w - H_w*(xmeasure - A_0*x0_prev - B_0*u0_prev);
        H_delta = -H_w*getD(x0_prev,u0_prev,A_1,A_2,B_1,B_2);
        % the Delta set of this time instant (H and h computed as in (3.3), (3.4))
        Delta_0 = Polyhedron(H_delta,h_delta);
        Delta_0.minHRep;
        %Delta_0.minVRep;
    end
    
    % Find $\Theta_t$
    Theta_t_M = Theta;
    for i = 0:(M-1)
        eval(['Set = Delta_' num2str(i) ';']);
        Theta_t_M = Theta_t_M & Set;
        Theta_t_M.minHRep;
        %Theta_t_M.minVRep;
    end
    %Theta_t = Theta_t_M & Theta_HC_prev;
    %H_thetat = Theta_t.A;
    %h_thetat = Theta_t.b;    
    H_thetat=[Theta_t_M.A;Theta_HC_prev.A];
    h_thetat=[Theta_t_M.b;Theta_HC_prev.b];
    Theta_t=Polyhedron(H_thetat,h_thetat);
    % Check if theta (true param vector) is in Theta_t
    check_belonging = H_thetat*theta_true - h_thetat;
    if check_belonging >0
        error('set estimation failed')
    end
    %% Find HC parameters 
    % find side
    for dim = 1 : p
        [~,lower(dim)]=linprog([zeros(dim-1,1);1;zeros(p-dim,1)],H_thetat,h_thetat,[],[],[],[],options);
        [~,upper(dim)]=linprog(-[zeros(dim-1,1);1;zeros(p-dim,1)],H_thetat,h_thetat,[],[],[],[],options);
        upper(dim) = -upper(dim);
    end
    HC_side = round(max(upper - lower),5);
    %projection
    HC_center = 0.5*(upper + lower)';
    %Projection_Set = theta_overline+0.5*(HC_side_prev-HC_side)*B_p;%Theta_HC_prev - 0.5*HC_side*B_p;
    for i=1:p       
       if HC_center(i)< HC_center_prev(i)-0.5*(HC_side_prev-HC_side)
           HC_center(i)=HC_center_prev(i)-0.5*(HC_side_prev-HC_side);
       elseif HC_center(i)>HC_center_prev(i)+0.5*(HC_side_prev-HC_side)
         HC_center(i)=HC_center_prev(i)+0.5*(HC_side_prev-HC_side);             
       else
         %no projection necessary  
       end          
    end
    % check that Theta^HC_t is in Theta^HC_t-1
    errorHC = norm(theta_true-HC_center,2);
    theta_HC_error = [theta_HC_error; errorHC];    
     Theta_HC=HC_center+HC_side*B_p;
   end
    %% Parameter estimation (MLS)
    if adapt==true
    if ii > 1
        theta_hat_prev = theta_RLS;
        [A_prev,B_prev] = getAB(theta_hat_prev,A_0,A_1,A_2,B_0,B_1,B_2);
        x0_est = A_prev*x0_prev + B_prev*u0_prev;
        theta_tilde_t = theta_hat_prev + mu_gain*getD(x0_prev,u0_prev,...
                                         A_1,A_2,B_1,B_2)'*(xmeasure-x0_est);
        [arg,~] = fmincon(@(x) norm(x-theta_tilde_t),zeros(p,1),Theta_HC.A,Theta_HC.b,...
                          [],[],[],[],[],options);%<- this comment can be implemented a lot more efficiently; just clipping
        theta_RLS = arg;
        %theta_RLS=theta_tilde_t;  
        if max(Theta_HC.A*theta_RLS-Theta_HC.b)>0
           error('LMS estimate not correctly projected') 
        end
        errorRLS = norm(theta_true-theta_RLS,2);
        theta_RLS_error = [theta_RLS_error; errorRLS];
        theta_RLS_rec=[theta_RLS_rec,theta_RLS];
    end
    end
   
    %% rho_t
   compute_rho(theta_overline,A_0,A_1,A_2,B_0,B_1,B_2,Hx,K);
    % L_Theta_t
    %% Start optimization
    
    theta_hat = theta_RLS;
    theta_overline = HC_center;
    eta = HC_side;
    
    % Set initial guess     
    %get init solution
    %reset beq with x_init
    b_eq(1:2*n)=[xmeasure;xmeasure];
    %compute A_eq,A_ineq using rho,eta
    A_eq_t=A_eq+...
        A_eq_theta_hat_1*theta_hat(1)+A_eq_theta_hat_2*theta_hat(2)+...
        A_eq_theta_overline_1*theta_overline(1)+A_eq_theta_overline_2*theta_overline(2);
    A_ineq_t=A_ineq+...
        A_ineq_rho*rho+...
        A_ineq_eta*eta;
    %update terminal set constraint with x_ss
    b_ineq_t=b_ineq+b_ineq_xss*x_ss;
    %update cost
    u_ss=[[0,1]-A_0(2,:)-A_1(2,:)*theta_hat(1)-A_2(2,:)*theta_hat(2)]*x_ss/B_0(2,1);
    y_s=[repmat(x_ss,2*(N+1),1);repmat(u_ss-K*x_ss,N,1);zeros(N+1,1)];
    %y_s=[repmat(x_ss,2*(N+1),1);repmat(u_ss,N,1);zeros(N+1,1)];
    f=-y_s'*H;
    t_Start = tic;
    [y_OL,V,exitflag] = quadprog(H,f,A_ineq_t,b_ineq_t,A_eq_t,b_eq,[],[],y_init,options);
    %<- more efficient implementation: use, e.g., acados and set online
    %changing variables as parameters
    t_Elapsed = toc(t_Start)
    if exitflag==-2
       error('problem infeasible') 
    end
    x_overline_OL = y_OL(1:n*(N+1));
    x_hat_OL=y_OL(n*(N+1)+1:2*n*(N+1));
    v_OL=y_OL(2*(N+1)*n+1:2*(N+1)*n+m*N);
    s_OL=y_OL(2*(N+1)*n+m*N+1:end);
    
    %% Store closed loop data
    t = [ t, tmeasure ];
    x = [ x, xmeasure ];
    x1 = [x1; xmeasure(1)];
    x1_ref = [x1_ref; x_ss(1)];
    u = [ u, v_OL(1:m)+K*xmeasure(1:n) ];
    s = [ s, s_OL(1) ];
    theta_LMS_store=[theta_LMS_store;theta_hat];
    theta_center_store=[theta_center_store;theta_overline];
    eta_store=[eta_store;eta];
    %% Save past state value for RLS
    x0_prev = xmeasure;
    u0_prev = u(:,end);
    
    %% Update closed-loop system (apply first control move to system)
    x_t = xmeasure;

    xmeasure = dynamic(x_t,u(:,end),theta_true,A_0,A_1,A_2,B_0,B_1,B_2,W_V);
    tmeasure = tmeasure + 1;
    if ii==1
       save('ARMPC_OL') 
    end
    %% Compute initial guess for next time step, based on terminal LQR controller (K_loc)
    
    %% Print and Plot
    % Plot predicted and closed-loop state trajetories
    f1 = figure(1);
    
    % State tube at time t = 0
    if mod(ii-1,step_width) == 0 
        %figure
        hold on
        plot(x_overline_OL(1:n:n*(N+1)),x_overline_OL(n:n:n*(N+1)),'g','Linewidth',1)
        plot(x_hat_OL(1:n:n*(N+1)),x_hat_OL(n:n:n*(N+1)),'b')
        %p1 = plot(x_OL(1:n:n*(N+1)),x_OL(n:n:n*(N+1)),'og');
        plot(x_ss(1),x_ss(2),'*')
        line([x_1max x_1max], [-3 3],'LineWidth',1.5)
        line([x_1min x_1min], [-3 3],'LineWidth',1.5)
        xlim([x_1min-0.2 x_1max+0.2])
        ylim([-3 3])
    elseif false%mod(ii-1,step_width) == step_width-1 
        %%
        figure(5)
        line([x_1max x_1max], [-3 3],'LineWidth',1.5)
        line([x_1min x_1min], [-3 3],'LineWidth',1.5)
        xlim([x_1min-0.2 x_1max+0.2])
        ylim([-3 3])
        plot(x_overline_OL(1:n:n*(N+1)),x_overline_OL(n:n:n*(N+1)),'g--*','Linewidth',1)
             hold on
        plot(x_hat_OL(1:n:n*(N+1)),x_hat_OL(n:n:n*(N+1)),'b--*')
        %%
        xmeasure
        %x_overline_OL(1:4)
        %x_hat_OL(1:4)
        %compute 
        %%
    end 
    % plot(x_OL(1:n:n*(N+1)),x_OL(n:n:n*(N+1)),'g')
    plot(x(1,:),x(2,:),'b'), grid on, hold on,
    p3 = plot(x(1,:),x(2,:),'xb');
    %legend([p1 p3],{'x OL t=0','x CL'})
    
    xlabel('x(1)')
    ylabel('x(2)')
    drawnow
end
MPCtime = toc(MPCstart);
hold off
%% PLOTS
figure(2)
stairs(t,u);
title('Input u')
figure(3)
plot(x1);
hold on
plot(x1_ref);
ylim([x_1min;x_1max])
title('Position x')
hold off
figure(4)
plot(theta_RLS_error);
title('RLS estimation error')
figure(5)
plot(theta_HC_error);
title('HC estimation error')
if adapt==true
save('ARMPC');
else
save('RMPC');
end
%% Display some values
disp('*****************************************************');
disp(['Total time for online computations: ',num2str(MPCtime)]);
disp(['Number of sides of X_0: ',num2str(c_0)]);
disp(['c_max = ',num2str(cmax)]);
disp(['d_bar = ',num2str(dbar)]);
sizey = size(y_init);
sizecon = size([A_ineq_t;A_eq_t]);
disp(['size(y) = ',num2str(sizey(1))]);
disp(['size(con) = ',num2str(sizecon(1))]);
if theta_out > 0
   error(['The true parameter was out of Theta_t ',num2str(theta_out),' times']);
end
end
%% FUNCTIONS
function xplus=dynamic(x,u,theta,A_0,A_1,A_2,B_0,B_1,B_2,W_V)
    [A,B] = getAB(theta,A_0,A_1,A_2,B_0,B_1,B_2);
    w = W_V(1,:)';
    xplus = A*x + B*u + (1-2*rand)*w;

end
function [A,B] = getAB(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function [A,B] = getABunc(p,A_1,A_2,B_1,B_2)
    A = p(1)*A_1 + p(2)*A_2;
    B = p(1)*B_1 + p(2)*B_2;
end

function [y] = getD(x,u,A_1,A_2,B_1,B_2)
    y = [A_1*x + B_1*u, A_2*x + B_2*u];
end 
 %%
 function rho=compute_rho(theta_overline,A_0,A_1,A_2,B_0,B_1,B_2,Hx,K) 
 options = optimset('Display','off');

    array_rho = [];
    [A_rho,B_rho] = getAB(theta_overline,A_0,A_1,A_2,B_0,B_1,B_2);
    A_K = A_rho + B_rho*K;
    c_0=size(Hx,1);
    for i = 1:c_0
        [~,fval] = linprog(-Hx(i,:)*A_K, Hx, ones(c_0,1),[],[],[],[],[],options);
        array_rho = [array_rho; -fval];
    end
    rho = max(array_rho);
 end