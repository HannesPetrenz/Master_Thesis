
function [P,K] = get_PK_LMIs(Q,R,Theta,rho_PK,A_0,A_1,A_2,B_0,B_1,B_2)
%LMIs to compte Lyapunpov P and feedback K for uncertain LTI system
% System dynamics
options = optimset('Display','off');

% System dimensions
n = length(Q(1,:)); % state dimension
m = length(R(1,:)); % input dimension

% Weighting matrices for cost function
Q_eps = sqrtm(Q);
R_eps = sqrt(R);

% Constraints on parameters and parameters variation
p1_min = min(Theta.V(:,1));
p1_max = max(Theta.V(:,1));
p2_min = min(Theta.V(:,2));
p2_max = max(Theta.V(:,2));

% X Y and Lambda matrices generation
Y_0_LMIs = sdpvar(m,n); 
X_0_LMIs = sdpvar(n);
% t = tic;

%% Multi Convexity constraints
% Starting to build the contraints
con = [];

% We build also constraints (21b) and (21c)
% The big one with matrices (2n+n+m)x(2n+n+m)
% And the one with Xmin<=X
i = 0;
% We check conditions only for the vertoces of ThetaBAR
for p1 = [p1_min,p1_max]
    for p2 = [p2_min,p2_max]

        % Parameter vector p
        p = [p1; p2];

        % Matrices X and Y
        Y = Y_0_LMIs; % + Y_1*p1 + Y_2*p2 + Y_3*p3;
        X = X_0_LMIs; % + X_1*p1 + X_2*p2 + X_3*p3;
        
        % We build matrices A and B dependent on p
        [A,B] = getAp(p,A_0,A_1,A_2,B_0,B_1,B_2);
        
        % We build the first matrix in (21b) and complete the constraints
        ineq = [X, (A*X+B*Y)', Q_eps*X, (R_eps*Y)';...
                A*X+B*Y, X, zeros(n,n+m);... 
                X*Q_eps, zeros(n), eye(n), zeros(n,m);...
                R_eps*Y, zeros(m,2*n), eye(m)];
        con = [con;ineq>=0];
        % Now 'con' is complete and we can start the optimization problem
        ineq = [rho_PK*X, (A*X+B*Y)';...
                A*X+B*Y, rho_PK*X];
        con = [con;ineq>=0];
        
        i=i+1;
    end
 end

% Display the number of iterations
disp('Number of iterations to get constraints for optimization to get P and K');
disp(i);
%% Optimization Problem
disp('Starting optimization');

% optimize(con,-log(det(X_max)))
optimize(con,-log(det(X_0_LMIs)));

Y_0_LMIs = value(Y_0_LMIs);
X_0_LMIs = value(X_0_LMIs);
% X_max = value(X_max);

P = inv(X_0_LMIs)
K = Y_0_LMIs*P
eigP = eig(inv(X_0_LMIs))

% toc(t)
% Store the output
save('convex_disc.mat','Y_0_LMIs','X_0_LMIs');


%% Check that P and K are fine
% ML matrices
% P = [1.467,0.207;0.207,1.731];
% K = [0.017,-0.41];

i=0;
for p1 = [p1_min,p1_max]
    for p2 = [p2_min,p2_max]
        % Parameter vector p
        p = [p1; p2; ];
        
        % We build matrices A and B dependent on p
        [A,B] = getAp(p,A_0,A_1,A_2,B_0,B_1,B_2);
        
        A_cl = A + B*K;
        
        M_check = A_cl.'*P*A_cl - P + Q + K.'*R*K;
%         M_check = A_cl.'*P*A_cl - P;
        E = eig(M_check);
        if round(E(1),6)>0 || round(E(2),6)>0
            i = i + 1;
        end
    end
end
 
if i>0
    disp('P and K are NOT fine!');
elseif i==0
    disp('P and K are fine');
end
%% FUNCTION
function [A,B]=getAp(p,A_0,A_1,A_2,B_0,B_1,B_2)
A = A_0 + p(1)*A_1 + p(2)*A_2 ;
B = B_0 + p(1)*B_1 + p(2)*B_2 ;
end
end