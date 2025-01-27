%This function plots the results from the main programm
close all
%clear
clc 
%Load the data
load("result_N8_iter10_deterministictrue.mat")

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
title("Inital solution and x^h")
p=[];
for h=1:length(X_RALMPC)
    p(h)=plot(X_RALMPC{h}(1,:),X_RALMPC{h}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
    hold on
    text{h}="X_{RALMPC,"+string(h)+"}";
end
legend(p,text)
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
for t=1:length(X_hat_OL_RALMPC{h})
    for k=1:length(X_hat_OL_RALMPC{h}{t})
        x_bar_OL=X_hat_OL_RALMPC{h}{t};
        s_OL=S_OL_RALMPC{h}{t};
        X_t=Polyhedron(H,s_OL(k)*ones(size(H,1),1)+H*x_bar_OL(:,k));
        plot(X_t,"color",[0,0,1])
        hold on
    end
    grid on
    p=plot(X_RALMPC{h}(1,:),X_RALMPC{h}(2,:),"Marker",'.', 'MarkerSize', 20);
end
legend(p,"closed loop solution")
%% Figure 3
figure(3)
Theta_HC0.plot('color', 'red');
hold on
Theta_HC_t_RALMPC.plot('color', 'lightblue');
xlabel("\theta_1")
ylabel("\theta_2")
title("Parameter set")
%% Figure 4
figure(4)
x=[];
y=[];
for i=1:length(J_RALMPC)
    x=[x,"RALMPC "+string(i)];
    y=[y,J_RALMPC{i}];
end
x= categorical(x);
bar(x,y)
grid on

