clear
clc
%search for all datasets
mat = dir('/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Data/*.mat'); 
%% loop through all data sets
Numberofiteration=[];
for q = 1:length(mat) 
    load(mat(q).folder+"/"+mat(q).name);
    plotten(H,X_bar_inital,S_inital,X_RALMPC,X_RMPC,X_RAMPC,X_OS,Theta_HC0,Theta_HC_t_RAMPC,Theta_HC_t_RALMPC,J_RMPC,J_RAMPC,J_RALMPC,J_RLMPC,J_OS,N,numberitertions,disturbance_deter,t_cpu_RMPC,t_cpu_RAMPC,t_cpu_RLMPC,t_cpu_RALMPC,t_cpu_OS,T_s)
    X_RALMPC_N{q}=X_RALMPC{end};
    J_RALMPC_N{q}=J_RALMPC{end};
    Numberofiteration(q)=N;
    T_cpu_RALMPC{q}=mean(vertcat(t_cpu_RALMPC{:}));
end
T_cpu_RAMPC=mean(vertcat(t_cpu_RAMPC{:}));
[id id] = sort(Numberofiteration,"descend");
Numberofiteration=Numberofiteration(id);
k=1;
for i=id
X_RALMPC_N_{k}=X_RALMPC_N{i};
J_RALMPC_N_{k}=J_RALMPC_N{i};
T_cpu_RALMPC_{k}=T_cpu_RALMPC{i};
k=k+1;
end
X_RALMPC_N=X_RALMPC_N_;
J_RALMPC_N=J_RALMPC_N_;
T_cpu_RALMPC=T_cpu_RALMPC_;
%% Additional Plots of the iterations
close all
%Figure 1
figure(1);
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
for h=1:length(X_RALMPC_N)
    p(h+1)=plot(X_RALMPC_N{h}(1,:),X_RALMPC_N{h}(2,:),"LineWidth",2,"Marker",'.', 'MarkerSize', 20);
    hold on
    text{h+1}="$x_{\mathrm{RALMPC},N="+string(Numberofiteration(h))+"}$";
end
p(h+2)=plot(X_OS(1,:),X_OS(2,:), 'LineWidth', 2, ...
    'Marker', '*', 'MarkerSize', 8, 'Color', 'b');
%description
grid on
xlabel("$x_{1}$",'Interpreter','latex')
ylabel("$x_{2}$",'Interpreter','latex')
text{h+2}="$x_{\mathrm{OS}}$";
legend(p,text,"Location",'eastoutside','Interpreter','latex')
%save
name1="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure11_"+"Iter"+string(numberitertions)+"Disturbance"+string(disturbance_deter)+".eps";
print('-depsc2', name1); % EPS format
%Figure 2
figure(2);
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
x1=["RAMPC"];
y1=[J_RAMPC{end}];
x2=[];
y2=[];
for h=1:length(J_RALMPC_N)
    x2=[x2,"RALMPC N="+string(Numberofiteration(h))];
    y2=[y2,J_RALMPC_N{h}];
end
x= categorical([x1,x2]);
x = reordercats(x,cellstr(x)');
bar(x,[y1,y2])
hold on
p=yline(J_OS,"LineWidth",2,"Color","red");
legend(p,"Optimal cost",'Interpreter','latex')
grid on
ylim([120,16])
set(gcf, 'Color', 'w');
name2="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure12_"+"Iter"+string(numberitertions)+"Disturbance"+string(disturbance_deter)+".eps";
print('-depsc2', name2); % EPS format
%% Figure 3
figure(3);
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
x1=["T_{RAMPC}"];
y1=[T_cpu_RAMPC];
x2=[];
y2=[];
for h=1:length(T_cpu_RALMPC)
    x2=[x2,"T_{RALMPC}(N="+string(Numberofiteration(h))+")"];
    y2=[y2,T_cpu_RALMPC{h}];
end
x= categorical([x1,x2]);
x = reordercats(x,cellstr(x)');
bar(x,[y1,y2])
grid on
title("CPU Times")
ylabel("CPU time [s]")
name3="/Users/hannes/Documents/MATLAB/Master_Thesis/RALMPC_LTI_Query/Scenario 2/Plots/Figure13_"+"Iter"+string(numberitertions)+"Disturbance"+string(disturbance_deter)+".eps";
exportgraphics(gcf, name3, 'ContentType', 'vector');

