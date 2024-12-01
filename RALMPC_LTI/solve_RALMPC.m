function [x_bar_OL,x_hat_OL,v_OL,s_OL,Q_wc]=solve_RALMPC(xmeasure,SS,Q_fun,H_cost,y_init,n,m,N,theta_bar_t,eta_t,rho_theta_t,theta_hat_t,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF,Q_wc,H,options)
J=inf;
x_bar_OL=[];
x_hat_OL=[];
v_OL=[];
s_OL=[];

%Decrease the computational burdens by decreasing the sample set;
idx=find(Q_fun<=Q_wc+5);
SS=SS(:,idx);
Q_fun=Q_fun(:,idx);
for j=1:size(SS,2)
    %Select xF and sF
    xF=SS(1:n,j);
    sF=SS(n+1,j);
    P=Q_fun(j);
    %Optimization
    %Update the constraints with the updated parameter
    [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(xmeasure,xF,sF,eta_t,rho_theta_t,theta_bar_t,theta_hat_t,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF,n,N);
    %Optimize the quadratic programm 
    [y_OL,V,exitflag] = quadprog(H_cost,zeros(n*(N+1)+n*(N+1)+m*N+(N+1)+1,1),A_ineqt,b_ineqt,A_eqt,b_eqt,[],[],y_init,options);
    if exitflag==1 && J>V+P
        %Extract the open loop solution
        x_bar_OL = reshape(y_OL(1:2*(N+1))',n,[]);
        x_hat_OL=reshape(y_OL(n*(N+1)+1:2*n*(N+1)),n,[]);
        v_OL=reshape(y_OL(2*(N+1)*n+1:2*(N+1)*n+m*N),m,[]);
        s_OL=y_OL(2*(N+1)*n+m*N+1:end-1);
        J=V+P;   
        Q_wc=P;
        %Update sample safe set
    end
end
if isempty(x_bar_OL)
    "error not feasible"
end
plotting(x_bar_OL,s_OL,H)
end


%% Help function
function [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(x0,xF,sF,eta,rho,theta_bar,theta_hat,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,b_ineq_xF,b_ineq_sF,n,N)
%This function builds the constraints for the current theta and rho.
    b_eq(1:n)=x0;
    b_eq(n*(N+1)+1:n*(N+2))=x0;
    A_eqt=A_eq+A_eq_theta1_bar*theta_bar(1)+A_eq_theta2_bar*theta_bar(2)+A_eq_theta1_hat*theta_hat(1)+A_eq_theta2_hat*theta_hat(2);
    b_eqt=b_eq;
    A_ineqt=A_ineq+A_ineq_eta*eta+A_ineq_rho*rho;
    b_ineqt=b_ineq+b_ineq_xF*xF+b_ineq_sF*sF;
end


function plotting(x_bar_OL,s_OL,H)
for k=1:length(x_bar_OL)       
        X_t=Polyhedron(H,s_OL(k)*ones(size(H,1),1)+H*x_bar_OL(:,k));
        plot(X_t,"color",[0 0 1])
        hold on
end
    plot(x_bar_OL(1,:),x_bar_OL(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20)
end
