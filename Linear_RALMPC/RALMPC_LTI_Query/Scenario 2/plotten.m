function plotten(H,X_bar_inital,S_inital,X_RALMPC,X_RMPC,X_RAMPC,X_OS,Theta_HC0,Theta_HC_t_RAMPC,Theta_HC_t_RALMPC,J_RMPC,J_RAMPC,J_RALMPC,J_RLMPC,J_OS,N,iteration,disturbance,t_cpu_RMPC,t_cpu_RAMPC,t_cpu_RLMPC,t_cpu_RALMPC,t_cpu_OS,Ts)
%This function plots the results from the main programm
close all
%% Figure 1
fig1 = figure(1);
%style
width = 16;     % Width in inches
height = 9;    % Height in inches
set(gcf,'InvertHardcopy','on');
pause(0.5);
set(gcf,'PaperUnits', 'inches');
pause(0.5);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(0.5);
%Plot the Polyhedron
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly,'Color', [0.8, 0.8, 1], 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    hold on
end
p(1)=plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{1}="$x_{\mathrm{intial}}$";
for h=1:length(X_RALMPC)
    p(h+1)=plot(X_RALMPC{h}(1,:),X_RALMPC{h}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
    hold on
    text{h+1}="$x_{\mathrm{RALMPC},"+string(h)+"}$";
end
p(h+2)=plot(X_RMPC(1,:),X_RMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","-.");
p(h+3)=plot(X_RAMPC{end}(1,:),X_RAMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","--");
p(h+4)=plot(X_OS(1,:),X_OS(2,:), 'LineWidth', 2, ...
    'Marker', '*', 'MarkerSize', 8, 'Color', 'b');
%description
grid on
xlabel("$x_{1}$",'Interpreter','latex')
ylabel("$x_{2}$",'Interpreter','latex')
title("Trajectory and intial solution",'Interpreter','latex')
text{h+2}="$x_{\mathrm{RMPC}}$";
text{h+3}="$x_{\mathrm{RAMPC}}$";
text{h+4}="$x_{\mathrm{OS}}$";
legend(p,text,"Location",'eastoutside','Interpreter','latex')
%save
name1="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure1_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".eps";
print('-depsc2', name1); % EPS format
%% Figure 2
figure(2);
%style
width = 16;     % Width in inches
height = 9;    % Height in inches
set(gcf,'InvertHardcopy','on');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%Plot the Polyhedron
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly,'Color', [0.8, 0.8, 1], 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    hold on
end
p=[];
text=[];
p(1)=plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{1}="$x_{\mathrm{intial}}$";
p(2)=plot(X_RALMPC{end}(1,:),X_RALMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{2}="$x_{\mathrm{RALMPC},"+string(h)+"}$";
p(3)=plot(X_RMPC(1,:),X_RMPC(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","-.");
p(4)=plot(X_RAMPC{end}(1,:),X_RAMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","--");
p(5)=plot(X_OS(1,:),X_OS(2,:), 'LineWidth', 2, ...
    'Marker', '*', 'MarkerSize', 8, 'Color', 'b');
%description
grid on
xlabel("$x_{1}$",'Interpreter','latex')
ylabel("$x_{2}$",'Interpreter','latex')
title("Trajectory and intial solution",'Interpreter','latex')
text{3}="$x_{\mathrm{RMPC}}$";
text{4}="$x_{\mathrm{RAMPC}}$";
text{5}="$x_{\mathrm{OS}}$";
legend(p,text,"Location",'eastoutside','Interpreter','latex')
%save
name2="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure2_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".eps";
print('-depsc2', name2); % EPS format
%% Figure 3
figure(3);
%style
width = 6;     % Width in inches
height =5;    % Height in inches
set(gcf,'InvertHardcopy','on');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
Theta_HC0.plot('color', 'red');
hold on
Theta_HC_t_RAMPC{end}.plot('color','lightgreen');
Theta_HC_t_RALMPC.plot('color', 'lightblue');
xlabel("$\theta_1$",'Interpreter','latex')
ylabel("$\theta_2$",'Interpreter','latex')
title("Parameter set",'Interpreter','latex')
legend("$\Theta^{\mathrm{HC}}_{0}$","$\Theta^{\mathrm{HC}}_{\mathrm{RAMPC}}$","$\Theta^{\mathrm{HC}}_{\mathrm{RALMPC}}$",'Interpreter','latex')
name3="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure3_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
print('-depsc2', name3); % EPS format
%% Figure 4
figure(4);
%style
width = 6;     % Width in inches
height =5;    % Height in inches
set(gcf, 'InvertHardCopy', 'off');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
x1=["RMPC","RAMPC"];
y1=[J_RMPC,J_RAMPC{end}];
x2=[];
y2=[];
x3=[];
y3=[];
for i=1:length(J_RALMPC)
    x3=[x3,"RALMPC "+string(i)];
    x2=[x2,"RLMPC "+string(i)];
    y3=[y3,J_RALMPC{i}];
    y2=[y2,J_RLMPC{i}];
end
x= categorical([x1,x2,x3]);
x = reordercats(x,cellstr(x)');
bar(x,[y1,y2,y3])
hold on
p=yline(J_OS,"LineWidth",2,"Color","red");
legend(p,"Optimal cost",'Interpreter','latex')
grid on
ylim([120,165])
set(gcf, 'Color', 'w');
name4="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure4_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
exportgraphics(gcf, name4, 'ContentType', 'vector');
%% Figure 5
figure(5);
%style
width = 10;     % Width in inches
height =8;    % Height in inches
set(gcf, 'InvertHardCopy', 'off');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
x=categorical(["RMPC","RAMPC","RLMPC","RALMPC"]);
x = reordercats(x,cellstr(x)');
y=[J_RMPC,J_RAMPC{end},J_RLMPC{end},J_RALMPC{end}];
bar(x,y)
hold on
p=yline(J_OS,"LineWidth",2,"Color","red");
legend(p,"Optimal Cost")
grid on
ylim([120,165])
title("Stationary Cost")
name5="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure5_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
exportgraphics(gcf, name5, 'ContentType', 'vector');
%% Figure 6
figure(6);
%style
width = 10;     % Width in inches
height =8;    % Height in inches
set(gcf, 'InvertHardCopy', 'off');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
x=[];
y=[];
for i=1:length(J_RALMPC)
    x=[x;i];
    y=[y;[J_RAMPC{i},J_RLMPC{i},J_RALMPC{i}]];
end
bar(x,y)
hold on
p=yline(J_OS,"LineWidth",2,"Color","red");
legend("RAMPC","RLMPC","RALMPC","Optimal Cost")
grid on
ylim([120,165])
title("Learning the Optimal Cost")
name6="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure6_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
exportgraphics(gcf, name6, 'ContentType', 'vector');
%% Figure 7
figure(7);
%style
width = 10;     % Width in inches
height =8;    % Height in inches
set(gcf, 'InvertHardCopy', 'off');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
T_cpu_RMPC=t_cpu_RMPC;
T_cpu_RAMPC=mean(vertcat(t_cpu_RAMPC{:}));
T_cpu_RLMPC=mean(vertcat(t_cpu_RLMPC{:}));
T_cpu_RALMPC=mean(vertcat(t_cpu_RALMPC{:}));
T_cpu_OS=t_cpu_OS;
x=categorical(["T_{RMPC}","T_{RAMPC}","T_{RLMPC}","T_{RALMPC}"]);
x = reordercats(x,cellstr(x)');
y=[T_cpu_RMPC,T_cpu_RAMPC,T_cpu_RLMPC,T_cpu_RALMPC];
bar(x,y)
grid on
title("CPU Times")
ylabel("CPU time [s]")
name7="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure7_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
exportgraphics(gcf, name7, 'ContentType', 'vector');
%% Figure 8
figure(8)
%style
width = 12;     % Width in inches
height =8;    % Height in inches
set(gcf, 'InvertHardCopy', 'off');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%plot
seq=0:length(X_OS(1,:))-1;
time=seq*Ts;
plot(time,X_OS(1,:),"LineWidth",2)
hold on
plot(time,X_RAMPC{end}(1,:),"LineWidth",2)
plot(time,X_RALMPC{end}(1,:),"LineWidth",2)
grid on
yline(0,"LineWidth",2,"Color","red");
legend("x_{1,OS}","x_{1,RAMPC}","x_{1,RALMPC}","Control Goal")
xlabel("time [s]")
ylabel("x_1")
name8="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure8_"+"Iter"+string(iteration)+"_N"+string(N)+"Determenistic"+string(disturbance)+".eps";
exportgraphics(gcf, name8, 'ContentType', 'vector');
%% Figure 9
figure(9);
%style
width = 16;     % Width in inches
height = 9;    % Height in inches
set(gcf,'InvertHardcopy','on');
pause(0.5);
set(gcf,'PaperUnits', 'inches');
pause(0.5);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(0.5);
%Plot the Polyhedron
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly,'Color', [0.8, 0.8, 1], 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    hold on
end
p=[];
text=[];
p(1)=plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{1}="$x_{\mathrm{intial}}$";
for h=1:length(X_RALMPC)
    p(h+1)=plot(X_RALMPC{h}(1,:),X_RALMPC{h}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
    hold on
    text{h+1}="$x_{\mathrm{RALMPC},"+string(h)+"}$";
end
p(h+2)=plot(X_RAMPC{end}(1,:),X_RAMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","--");
p(h+3)=plot(X_OS(1,:),X_OS(2,:), 'LineWidth', 2, ...
    'Marker', '*', 'MarkerSize', 8, 'Color', 'b');
%description
grid on
xlabel("$x_{1}$",'Interpreter','latex')
ylabel("$x_{2}$",'Interpreter','latex')
title("Trajectory and intial solution",'Interpreter','latex')
text{h+2}="$x_{\mathrm{RAMPC}}$";
text{h+3}="$x_{\mathrm{OS}}$";
legend(p,text,"Location",'eastoutside','Interpreter','latex')
%save
name9="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure9_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".eps";
print('-depsc2', name9); % EPS format
%% Figure 10
figure(10);
%style
width = 16;     % Width in inches
height = 9;    % Height in inches
set(gcf,'InvertHardcopy','on');
pause(1);
set(gcf,'PaperUnits', 'inches');
pause(1);
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
pause(1);
%Plot the Polyhedron
for k=1:length(X_bar_inital)
    X_poly=Polyhedron(H,S_inital(k)*ones(size(H,1),1)+H*X_bar_inital(:,k));
    plot(X_poly,'Color', [0.8, 0.8, 1], 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    hold on
end
p=[];
text=[];
p(1)=plot(X_bar_inital(1,:),X_bar_inital(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{1}="$x_{\mathrm{intial}}$";
p(2)=plot(X_RALMPC{end}(1,:),X_RALMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
text{2}="$x_{\mathrm{RALMPC},"+string(h)+"}$";
p(3)=plot(X_RAMPC{end}(1,:),X_RAMPC{end}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20,"LineStyle","--");
p(4)=plot(X_OS(1,:),X_OS(2,:), 'LineWidth', 2, ...
    'Marker', '*', 'MarkerSize', 8, 'Color', 'b');
%description
grid on
xlabel("$x_{1}$",'Interpreter','latex')
ylabel("$x_{2}$",'Interpreter','latex')
title("Trajectory and intial solution",'Interpreter','latex')
text{3}="$x_{\mathrm{RAMPC}}$";
text{4}="$x_{\mathrm{OS}}$";
legend(p,text,"Location",'eastoutside','Interpreter','latex')
%save
name10="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure10_"+"Iter"+string(iteration)+"_N"+string(N)+"Disturbance"+string(disturbance)+".eps";
print('-depsc2', name10); % EPS format
end