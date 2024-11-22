function [X_0] = lambda_contractive(lambda,F,G,A_K1,A_K2,A_K3,A_K4,K)
%compute lambda contractive set
options = optimset('Display','off');

L = F + G*K;    % The constrait set matrix for the state, in closed loop
V = L;
p = size(L,1);  % should be the total number of constraints we have (6 here)

nu = 0;
A_K = [];
fmax = -inf;
notFinished = true;
while(notFinished)
    for j = 1:4
        eval(['A_K=A_K' num2str(j) '/lambda;']);
        if nu == 0
            for i = 1:p
                [~,fval] = linprog(-V(end-p+i,:)*A_K, V, ones(size(V,1),1),[],[],[],[],[],options);
                fmax = max(-fval-1, fmax);
            end
        elseif nu > 0
            for i = 1:p*4
                [~,fval] = linprog(-V(end-p*4+i,:)*A_K, V, ones(size(V,1),1),[],[],[],[],[],options);
                fmax = max(-fval-1, fmax);
            end    
        end
    end
    
    if (fmax <= 0)           
        notFinished = 0;
    else
        fmax = -inf;
        Vend = [];
        for j = 1:4
            eval(['A_K=A_K' num2str(j) '/lambda;']);
            Vend = [Vend; V(end-p+1:end,:)*A_K];              % To V we add A_k^i*V (another "constraint")
        end
        V = [V; Vend];
    end
    nu=nu+1;         % Number of iterations to find the MRPI set
    disp(nu);
end
disp('Number of iterations to find X_0:');
disp(nu);
X_0 = Polyhedron(V,ones(length(V(:,1)),1));
end