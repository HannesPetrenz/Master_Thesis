function [SS_0,J_wc_0,X_bar,S,J_0]=computeinitialsolution(x0,x_s,u_s,theta_bar_0,rho_theta_0,eta_0,L_B_rho,d_bar,L,l,c,c_xs,P,mpciterations,Q,R,m,c_alpha)
X_bar=[x0];
S=[0];
U_bar=[];
u_bar_inital=[-2.0000,-1.9968,-1.8383,-1.1244,-0.4645,-0.0853,-0.1017,-0.1192,-0.1380,0.3518,1.0436,1.7076];
k=1;
%Compute inital solution 
while k<mpciterations
    if k<=length(u_bar_inital)
        u_bar_k=u_bar_inital(k);
    else 
        K_s = controller_K(u_s,x_s(1),x_s(2));
        u_bar_k=K_s*(X_bar(end)-x_s)+u_s;
    end
    U_bar=[U_bar,u_bar_k];
    [x_bar_kplus,s_kplus]=dynamics(X_bar(:,end),S(:,end),U_bar(:,end),theta_bar_0,rho_theta_0,eta_0,L_B_rho,d_bar);
    X_bar=[X_bar,x_bar_kplus];
    S=[S,s_kplus];
    k=k+1;
end
u_bar_N=K_s*(X_bar(end)-x_s)+u_s;
U_bar=[U_bar,u_bar_N];
%check constrain satisfaction 
issatisfied=true;
for i=1:size(U_bar,2)
    x_bar_k=X_bar(:,i);
    u_bar_k=U_bar(:,i);
    s_k=S(:,i);
    if ~all(L*[u_bar_k;x_bar_k]-l+c*s_k<=1e6)
        issatisfied=false;
    end
end
x_bar_N=X_bar(:,end);
s_N=S(:,end);
if ~(sqrt((x_bar_N-x_s).'*P*(x_bar_N-x_s))+s_N<=c_xs)
    issatisfied=false;
end    
if ~issatisfied
    fprintf("\tThe inital condition does not satisfy the constraints\n");
end
%Compute worst coste to go
cost_to_go_wc=c_alpha*worstcost(X_bar(:,end),0,P,S(:,end),P,0,[0,0],L(:,m+1:end),l); %terminal worst case cost
for i=size(X_bar,2)-1:-1:2
    K = controller_K(U_bar(:,i),X_bar(1,i),X_bar(2,i));
    cost=worstcost(X_bar(:,i),U_bar(:,i),P,S(:,i),Q,R,K,L(:,m+1:end),l);
    cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
end
J_wc_0=flip(cost_to_go_wc);
%Initilize the sample set
SS_0=[X_bar(:,2:end);S(:,2:end);U_bar(:,2:end)];
%Compute reference cost
J_0=0;
for i=1:length(X_bar)
    u_hat_k=U_bar(i);
    J_0=J_0+X_bar(:,i)'*Q*X_bar(:,i)+u_hat_k'*R*u_hat_k;
end
end

function [x_kplus_bar,s_kplus]=dynamics(x_bar_k,s_k,u_bar,theta_bar_0,rho_theta_0,eta_0,L_B_rho,d_bar)
    %system dynamic for theta_bar
    x_kplus_bar = system_f(theta_bar_0(1),theta_bar_0(2),u_bar,x_bar_k(1),x_bar_k(2));
    %tube dynamic
    w_Theta = max(uncertainty_w_Theta(eta_0,x_bar_k(1),x_bar_k(2)));
    w_delta_Theta_D=w_Theta+(eta_0*L_B_rho)*s_k+d_bar;
    s_kplus=rho_theta_0*s_k+w_delta_Theta_D;
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
