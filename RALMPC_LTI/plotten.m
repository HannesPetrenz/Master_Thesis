%This function plots the results from the main programm

%Load the data
load("result.mat")

%% Figure 1
figure(1)
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly)
    hold on
end
plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20)
grid on
xlabel("x_{1}")
ylabel("x_{2}")
title("Inital solution and x^h_t")
p=[];
for h=1:length(X)
    p(h)=plot(X{h}(1,:),X{h}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
    hold on   
end
p(h+1)=plot(X_RMPC(1,:),X_RMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
p(h+2)=plot(X_RAMPC(1,:),X_RAMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
p(h+3)=plot(X_OP(1,:),X_OP(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);

legend(p)
%% Figure 2
figure(2)
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly)
    hold on
end
plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20)
grid on
xlabel("x_{1}")
ylabel("x_{2}")
title("Inital solution and X_{k|t}^h")
h=3;
for t=1:length(X_hat_OL{h})
    for k=1:length(X_hat_OL{h}{t})
        x_bar_OL=X_hat_OL{h}{t};
        s_OL=S_OL{h}{t};
        X_t=Polyhedron(H,s_OL(k)*ones(size(H,1),1)+H*x_bar_OL(:,k));
        plot(X_t)
        hold on
    end
    grid on
    plot(X{h}(1,:),X{h}(2,:),"Marker",'.', 'MarkerSize', 20)
end
%% Figure 3
figure(3)
Theta_HC0.plot('color', 'red');
hold on
Theta_HC_t.plot('color', 'lightblue');
%% Figure 4
figure(4)
bar(1,J_RMPC,"LineWidth",2)
bar(2,J_RAMPC,"LineWidth",2)
 for i=1:length(J)
    bar(i+2,J{i})
    hold on
 end
 yline(J_OP,"LineWidth",2,"Color","red")

