function plotten(H,X_bar_inital,S_inital,X_RALMPC,X_RMPC,X_RAMPC,X_OS,Theta_HC0,Theta_HC_t_RALMPC,J_RMPC,J_RAMPC,J_RALMPC,J_RLMPC,J_OS,N,iteration,disturbance)
%This function plots the results from the main programm
close all
%% Figure 1
fig1 = figure(1);
fig1.WindowState = 'maximized';
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
p(h+1)=plot(X_RMPC(1,:),X_RMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","-.");
p(h+2)=plot(X_RAMPC(1,:),X_RAMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","-.");
p(h+3)=plot(X_OS(1,:),X_OS(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","-.");
text{h+1}="X_{RMPC}";
text{h+2}="X_{RAMPC}";
text{h+3}="X_{OS}";
legend(p,text,"Location",'eastoutside')
name1="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure1_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".png";
saveas(gcf,name1)
% %% Figure 2
% figure(2)
% for k=1:length(X_bar_inital)
%     X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
%     plot(X_poly)
%     hold on
% end
% plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20)
% grid on
% xlabel("x_{1}")
% ylabel("x_{2}")
% title("Inital solution and X_{k|t}^h")
% h=3;
% for t=1:length(X_hat_OL_RALMPC{h})
%     for k=1:length(X_hat_OL_RALMPC{h}{t})
%         x_bar_OL=X_hat_OL_RALMPC{h}{t};
%         s_OL=S_OL_RALMPC{h}{t};
%         X_t=Polyhedron(H,s_OL(k)*ones(size(H,1),1)+H*x_bar_OL(:,k));
%         plot(X_t,"color",[0,0,1])
%         hold on
%     end
%     grid on
%     p=plot(X_RALMPC{h}(1,:),X_RALMPC{h}(2,:),"Marker",'.', 'MarkerSize', 20);
% end
% legend(p,"closed loop solution")
%% Figure 3
figure(3);
Theta_HC0.plot('color', 'red');
hold on
Theta_HC_t_RALMPC.plot('color', 'lightblue');
xlabel("\theta_1")
ylabel("\theta_2")
title("Parameter set")
name2="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure2_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".png";
saveas(gcf,name2)
%% Figure 4
figure(4);
x1=["RMPC","RAMPC"];
y1=[J_RMPC,J_RAMPC];
x2=[];
y2=[];
x3=[];
y3=[];
for i=1:length(J_RALMPC)
    x2=[x2,"RALMPC "+string(i)];
    x3=[x3,"RLMPC "+string(i)];
    y2=[y2,J_RALMPC{i}];
    y3=[y3,J_RLMPC{i}];
end
x= categorical([x1,x2,x3]);
x = reordercats(x,cellstr(x)');
bar(x,[y1,y2,y3])
hold on
p=yline(J_OS,"LineWidth",2,"Color","red");
legend(p,"Optimal cost")
grid on
ylim([120,160])
name3="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure3_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".png";
saveas(gcf,name3)
end