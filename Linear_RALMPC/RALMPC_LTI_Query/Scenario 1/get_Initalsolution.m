function [SS_0,J_wc_0,X_bar,S]=get_Initalsolution(x0,s0,theta_bar0,B_p,eta_0,rho_theta0,L_B,d_bar,c,c_max,H,A_0,A_1,A_2,B_0,B_1,B_2,K,Q,R,P,F,G,m,n,p)
%compute the system matrix for theta_0
[A_cl_thetabar0,B_thetabar0] = get_AB_theta(theta_bar0,A_0,A_1,A_2,B_0,B_1,B_2);
A_cl_thetahat0=A_cl_thetabar0;
B_thetahat0=B_thetabar0;
%Simulate the dynamics
X_bar=[x0];
X_hat=[x0];
S=[s0];
V=[11	8	7	5	3	2	1	-0.5 -0.9	-1.4 -1.2 -1 -0 1.2];
k=1;
while k<20
    if k<=length(V)
        v=V(k);
    else 
        v=0;
    end
    [x_kplus_bar,x_kplus_hat,s_kplus]=dynamics(X_bar(:,k),X_hat(:,k),S(:,k),v,A_cl_thetabar0,B_thetabar0,A_cl_thetahat0,B_thetahat0,K,B_p,eta_0,rho_theta0,L_B,d_bar,H,A_1,A_2,B_1,B_2,p);
    X_bar=[X_bar,x_kplus_bar];
    X_hat=[X_hat,x_kplus_hat];
    S=[S,s_kplus];
    terminate=check_terminalcondition(X_bar(:,k),S(:,k),c_max,H);
    k=k+1;
end

%Compute the worst-case cost-to-go
%terminal cost
X=Polyhedron(H,S(end)*ones(size(H,1),1)+H*X_bar(:,end));
cost=0;
for i=1:size(H,1)
    cost=max([cost;X.V(i,:)*P*X.V(i,:)']);
end
cost_to_go_wc=cost;
%stage cost
for k=length(X_bar)-1:-1:2
    X=Polyhedron(H,S(k)*ones(size(H,1),1)+H*X_bar(:,k));
    cost=0;
    for i=1:size(H,1)
    cost=max([cost;X.V(i,:)*Q*X.V(i,:)'+(X.V(i,:))*K'*R*K*X.V(i,:)']);
    end
    cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
end
%Check constraints
issatisfied=true;
v=[];
for k=1:length(X_bar)
    if k<=length(V)
        v(k)=V(k);
    else 
        v(k)=0;
    end
    
    if any((F+G*K)*X_bar(:,k)+G*v+c*S(k)-1>0)
        issatisfied=false;
    end
end
if issatisfied==false
    disp("The inital solution violates the robust constraints")
end
%Construct inital sample set: 
J_wc_0=flip(cost_to_go_wc);
SS_0=[X_bar(:,2:end);S(:,2:end);v(:,2:end)];

checksecondcost=0;
for i=1:length(X_bar)
    u_second=v(i)+K*X_bar(:,i);
    checksecondcost=checksecondcost+X_bar(:,i)'*Q*X_bar(:,i)+u_second'*R*u_second;
end
end


%% Help function
function [A,B] = get_AB_theta(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function [x_kplus_bar,x_kplus_hat,s_kplus]=dynamics(x_k_bar,x_k_hat,s_k,v_k,A_cl_thetabar0,B_thetabar0,A_cl_thetahat0,B_thetahat0,K,B_p,eta_0,rho_theta0,L_B,d_bar,H,A_1,A_2,B_1,B_2,p)
%Input
u_k_hat=v_k+K*x_k_hat;
u_k_bar=v_k+K*x_k_bar;
%system dynamic
x_kplus_bar=A_cl_thetabar0*x_k_bar+B_thetabar0*u_k_bar;
x_kplus_hat=A_cl_thetahat0*x_k_hat+B_thetahat0*u_k_hat;
%Overapproximate the uncertainty
w_k=get_overappuncertainty(x_k_bar,u_k_bar,s_k,B_p,eta_0,rho_theta0,L_B,H,d_bar,A_1,A_2,B_1,B_2,p);
%scalar tube dynamics
s_kplus=rho_theta0*s_k+w_k;
end

function w_k=get_overappuncertainty(x_k_bar,u_k_bar,s_k,B_p,eta_0,rho_theta0,L_B,H,d_bar,A_1,A_2,B_1,B_2,p)
w_array=[];
D=get_D(x_k_bar,u_k_bar,A_1,A_2,B_1,B_2);
for l=1:2^p
    for i=1:size(H,1)
    w_array=[w_array;d_bar+eta_0*(L_B*s_k+H(i,:)*D*B_p.V(l,:)')];
    end
end
w_k=max(w_array);
end

function D=get_D(x_k_bar,u_k_bar,A_1,A_2,B_1,B_2)
D=[A_1*x_k_bar+B_1*u_k_bar,A_2*x_k_bar+B_2*u_k_bar];
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