clear
clc
close all
%% Define symbolic variables
syms m Iz lf lr mu g
syms v_x v_y psi psi_dot s e_y e_psi e_psi_dot
syms s_k e_y_k e_psi_k r_k v_x_k v_y_k
syms delta_f a
syms delta_f_k a_k T
syms C D B
syms r
%% Load track data
track=load("L_track_barc.mat");
[x_grid,y_grid,track_key_pts]=racetrack(track.cl_segs, [0,0,0],track.slack,track.track_width);
plot(x_grid(:,1),y_grid(:,1))
hold on
plot(x_grid(:,2),y_grid(:,2))
%% Tire forces
%slip angles
alpha_f=atan((v_y+lf*psi_dot)/norm(v_x));
alpha_r=atan((v_y-lr*psi_dot)/norm(v_x));
alpha_f_k=atan((v_y_k+lf*r_k)/norm(v_x_k));
alpha_r_k=atan((v_y_k-lr*r_k)/norm(v_x_k));
%Parcejka formula
fp_f=D*sin(C*atan(B*alpha_f));
fp_r=D*sin(C*atan(B*alpha_r));
fp_f_k=D*sin(C*atan(B*alpha_f_k));
fp_r_k=D*sin(C*atan(B*alpha_r_k));
%force
F_yf=-0.5*m*g*mu*fp_f;
F_yr=-0.5*m*g*mu*fp_r;
F_yf_k=-0.5*m*g*mu*fp_f_k;
F_yr_k=-0.5*m*g*mu*fp_r_k;
%% Curvature of the centerline
kappa=1/r;
%% Dynamic equation
s_dot=(v_x*cos(e_psi)-v_y*sin(e_psi))/(1-e_y*kappa);
e_y_dot=v_x*sin(e_psi)+v_y*cos(e_psi);
e_psi_dot=psi_dot-kappa*(v_x*cos(e_psi)-v_y*sin(e_psi))/(1-e_y*kappa);
psi_ddot=1/Iz*(lf*F_yf-lr*F_yr);
v_x_dot=a+psi_dot*v_y;
v_y_dot=1/m*(F_yf*cos(delta_f)+F_yr)-psi_dot*v_x;
%system
f=[s_dot;e_y_dot;e_psi;psi_ddot;v_x_dot;v_y_dot];
x=[s;e_y;e_psi;psi_dot;v_x;v_y];
u=[a;delta_f];
%% Discretization
s_kplus=s_k+T*(v_x_k*cos(e_psi_k)-v_y_k*sin(e_psi_k))/(1-e_y_k*kappa);
e_y_kplus=e_y_k+T*(v_x_k*sin(e_psi_k)+v_y_k*cos(e_psi_k));
e_psi_kplus=e_psi_k+T*(r_k-kappa*(v_x_k*cos(e_psi_k)-v_y_k*sin(e_psi_k))/(1-e_y_k*kappa));
r_kplus=r_k+T*(1/Iz*(lf*F_yf_k-lr*F_yr_k));
v_x_kplus=v_x_k+T*(a_k+r_k*v_y_k);
v_y_kplus=v_y_k+T*(1/m*(F_yf_k*cos(delta_f_k)+F_yr_k)-psi_dot*v_x_k);
%system
f_k=[s_kplus;e_y_kplus;e_psi_kplus;r_kplus;v_x_kplus;v_y_kplus];
x_k=[s_k;e_y_k;e_psi_k;r_k;v_x_k;v_y_k];
u_k=[a_k,delta_f_k];
%% Linearization
A_matrix=simplify(jacobian(f,x));
B_matrix=simplify(jacobian(f,u));
A_k=simplify(jacobian(f_k,x_k));
B_k=simplify(jacobian(f_k,u_k));
%% Sensitivity analysis 
s_rang=[min(track_key_pts(:,4)),max(track_key_pts(:,4))];
e_y_rang=[-track.track_width/2,track.track_width/2];
e_psi_rang=[-65/180*pi,65/180*pi];
psi_dot_rang=[0,10];
v_x_rang=[0.1,3];
v_y_rang=[-0.5,0.5];
delta_f_rang=[-0.4,0.4];%rad
a_rang=[-1.3,3];%m/s^2
m_=1.75;%kg
lf_=0.125;%m
lr_=0.125;%m
Iz_=0.03;%kgm^2
g_=9.81;%m/s^2
mu_=0.85;
B_=6;
C_=1.6;
D_=1;
T=0.1;
A_3D=[];
B_3D=[];
i=1;
for delta_f_=delta_f_rang
    for r_=unique(track.cl_segs(:,2)')
        for s_=linspace(s_rang(1),s_rang(2),4)
            for e_y_=linspace(e_y_rang(1),e_y_rang(2),4)
                for e_psi_=linspace(e_psi_rang(1),e_psi_rang(2),4)
                    for psi_dot_=linspace(psi_dot_rang(1),psi_dot_rang(2),4)
                        for v_x_=linspace(v_x_rang(1),v_x_rang(2),4)
                            for v_y_=linspace(v_y_rang(1),v_y_rang(2),4)
                                A_3D(:,:,i)=vpa(subs(A_matrix,[m,lf,lr,Iz,g,mu,B,C,D,s,e_y,e_psi,psi_dot,v_x,v_y,r,delta_f],[m_,lf_,lr_,Iz_,g_,mu_,B_,C_,D_,s_,e_y_,e_psi_,psi_dot_,v_x_,v_y_,r_,delta_f_]),3);
                                B_3D(:,:,i)=vpa(subs(B_matrix,[m,lf,lr,Iz,g,mu,B,C,D,s,e_y,e_psi,psi_dot,v_x,v_y,r,delta_f],[m_,lf_,lr_,Iz_,g_,mu_,B_,C_,D_,s_,e_y_,e_psi_,psi_dot_,v_x_,v_y_,r_,delta_f_]),3);                         
                                i=i+1;
                            end
                        end
                    end    
                end
            end
        end
    end
end
%%
A_max=[max(A_3D(1,1,:)),max(A_3D(1,2,:)),max(A_3D(1,3,:)),max(A_3D(1,4,:)),max(A_3D(1,5,:)),max(A_3D(1,6,:));
    max(A_3D(2,1,:)),max(A_3D(2,2,:)),max(A_3D(2,3,:)),max(A_3D(2,4,:)),max(A_3D(2,5,:)),max(A_3D(2,6,:));
    max(A_3D(3,1,:)),max(A_3D(3,2,:)),max(A_3D(3,3,:)),max(A_3D(3,4,:)),max(A_3D(3,5,:)),max(A_3D(3,6,:));
    max(A_3D(4,1,:)),max(A_3D(4,2,:)),max(A_3D(4,3,:)),max(A_3D(4,4,:)),max(A_3D(4,5,:)),max(A_3D(4,6,:));
    max(A_3D(5,1,:)),max(A_3D(5,2,:)),max(A_3D(5,3,:)),max(A_3D(5,4,:)),max(A_3D(5,5,:)),max(A_3D(5,6,:));
    max(A_3D(6,1,:)),max(A_3D(6,2,:)),max(A_3D(6,3,:)),max(A_3D(6,4,:)),max(A_3D(6,5,:)),max(A_3D(6,6,:))]
A_min=[min(A_3D(1,1,:)),min(A_3D(1,2,:)),min(A_3D(1,3,:)),min(A_3D(1,4,:)),min(A_3D(1,5,:)),min(A_3D(1,6,:));
    min(A_3D(2,1,:)),min(A_3D(2,2,:)),min(A_3D(2,3,:)),min(A_3D(2,4,:)),min(A_3D(2,5,:)),min(A_3D(2,6,:));
    min(A_3D(3,1,:)),min(A_3D(3,2,:)),min(A_3D(3,3,:)),min(A_3D(3,4,:)),min(A_3D(3,5,:)),min(A_3D(3,6,:));
    min(A_3D(4,1,:)),min(A_3D(4,2,:)),min(A_3D(4,3,:)),min(A_3D(4,4,:)),min(A_3D(4,5,:)),min(A_3D(4,6,:));
    min(A_3D(5,1,:)),min(A_3D(5,2,:)),min(A_3D(5,3,:)),min(A_3D(5,4,:)),min(A_3D(5,5,:)),min(A_3D(5,6,:));
    min(A_3D(6,1,:)),min(A_3D(6,2,:)),min(A_3D(6,3,:)),min(A_3D(6,4,:)),min(A_3D(6,5,:)),min(A_3D(6,6,:))]
B_max=[max(B_3D(1,1,:)),max(B_3D(1,2,:));
    max(B_3D(2,1,:)),max(B_3D(2,2,:));
    max(B_3D(3,1,:)),max(B_3D(3,2,:));
    max(B_3D(4,1,:)),max(B_3D(4,2,:));
    max(B_3D(5,1,:)),max(B_3D(5,2,:));
    max(B_3D(6,1,:)),max(B_3D(6,2,:))]
B_min=[min(B_3D(1,1,:)),min(B_3D(1,2,:));
    min(B_3D(2,1,:)),min(B_3D(2,2,:));
    min(B_3D(3,1,:)),min(B_3D(3,2,:));
    min(B_3D(4,1,:)),min(B_3D(4,2,:));
    min(B_3D(5,1,:)),min(B_3D(5,2,:));
    min(B_3D(6,1,:)),min(B_3D(6,2,:))]

%% Save
save('Analyse_Matrix.mat',"A_max","A_min","B_max","B_min")
%% Help function
function [x_grid,y_grid,track_key_pts]=racetrack(cl_segs, init_pos,slack,track_width)
half_width=track_width/2;
track_key_pts=get_track_key_pts(cl_segs, init_pos);
track_length = track_key_pts(end, 4);
seg_x = track_key_pts(:, 1);
seg_y = track_key_pts(:, 2);
seg_t = track_key_pts(:, 3);
cum_l = track_key_pts(:, 4);
seg_l = track_key_pts(:, 5);
seg_c = track_key_pts(:, 6);
%Cumulative track heading
cum_t = seg_t(1);
for i = 1:length(track_key_pts)-1 
    cum_t=[cum_t;cum_t(end) + seg_l(i+1)*seg_c(i+1)];
end
%Cumulative change in track heading
cum_dt = cum_t - cum_t(1);
%Get the x-y extents of the track
s_grid = linspace(0, track_length, 10 * track_length);
x_grid=[]; 
y_grid = [];
for i=1:length(s_grid)
    [xp, yp, psi] = local_to_global([s_grid(i), half_width + slack, 0],track_length,track_key_pts);
    [xm, ym, psi] = local_to_global([s_grid(i), -half_width - slack, 0],track_length,track_key_pts);
    x_grid=[x_grid;[xp,xm]];
    y_grid=[y_grid;[yp,ym]];
end
end

function track_key_pts=get_track_key_pts(cl_segs, init_pos)
n_segs = length(cl_segs);

track_key_pts = zeros(n_segs + 1, 6);
track_key_pts(1, 1) = init_pos(1);
track_key_pts(1, 2) = init_pos(2);
track_key_pts(1, 3) = init_pos(3);
for i=2:n_segs + 1
    x_prev = track_key_pts(i - 1, 1);
    y_prev = track_key_pts(i - 1, 2);
    psi_prev = track_key_pts(i - 1, 3);
    cum_s_prev = track_key_pts(i - 1, 4);

    l = cl_segs(i - 1, 1);
    r = cl_segs(i - 1, 2);

    if r == 0
        %No curvature (i.e. straight line)
        psi = psi_prev;
        x = x_prev + l * cos(psi_prev);
        y = y_prev + l * sin(psi_prev);
        curvature = 0;
    else
        %Find coordinates for center of curved segment
        x_c = x_prev - r * (sin(psi_prev));
        y_c = y_prev + r * (cos(psi_prev));
        %Angle spanned by curved segment
        theta = l / r;
        %end of curve
        x = x_c + r * sin(psi_prev + theta);
        y = y_c - r * cos(psi_prev + theta);
        %curvature of segment
        curvature = 1 / r;

        %next angle
        psi = wrap_angle(psi_prev + theta);
    end
    cum_s = cum_s_prev + l;
    track_key_pts(i,:) = [x, y, psi, cum_s, l, curvature];
end    
end

function wrapped_angle=wrap_angle(theta)
    if theta < -pi
        wrapped_angle = 2 * pi + theta;
    elseif theta > pi
        wrapped_angle = theta - 2 * pi;
    else
        wrapped_angle = theta;
    end
end


function [x, y, psi]=local_to_global(cl_coord,track_length,key_pts)
    s = cl_coord(1);
    while s < 0  
        s =s+track_length;
    end
    while s >= track_length
        s = s-track_length;
    end
    
    e_y = cl_coord(2);
    e_psi = cl_coord(3);
    
    %Find key point indicies corresponding to current segment
    %key_pts = [x, y, psi, cumulative length, segment length, signed curvature]
    key_pt_idx_s = find(s >= key_pts(:, 4));
    key_pt_idx_s=key_pt_idx_s(end);
    key_pt_idx_f = key_pt_idx_s + 1;
    seg_idx = key_pt_idx_s;
    
    x_s = key_pts(key_pt_idx_s, 1);
    y_s = key_pts(key_pt_idx_s, 2);
    psi_s = key_pts(key_pt_idx_s, 3);
    curve_s = key_pts(key_pt_idx_s, 6);
    x_f = key_pts(key_pt_idx_f, 1);
    y_f = key_pts(key_pt_idx_f, 2);
    psi_f = key_pts(key_pt_idx_f, 3);
    curve_f = key_pts(key_pt_idx_f, 6);
    
    l = key_pts(key_pt_idx_f, 5);
    d = s - key_pts(key_pt_idx_s, 4);  %Distance along current segment
    
    if curve_f == 0
        %Segment is a straight line
        x = x_s + (x_f - x_s) * d / l + e_y * cos(psi_f + pi / 2);
        y = y_s + (y_f - y_s) * d / l + e_y * sin(psi_f + pi / 2);
        psi = wrap_angle(psi_f + e_psi);
    else
        r = 1 / curve_f;
        dir = sign(r);
    
        %Find coordinates for center of curved segment
        x_c = x_s + abs(r) * cos(psi_s + dir * pi / 2);
        y_c = y_s + abs(r) * sin(psi_s + dir * pi / 2);
    
        %Angle spanned up to current location along segment
        span_ang = d / abs(r);
    
        %Angle of the tangent vector at the current location
        psi_d = wrap_angle(psi_s + dir * span_ang);
    
        ang_norm = wrap_angle(psi_s + dir * pi / 2);
        ang = -sign(ang_norm) * (pi - abs(ang_norm));
    
        x = x_c + (abs(r) - dir * e_y) * cos(ang + dir * span_ang);
        y = y_c + (abs(r) - dir * e_y) * sin(ang + dir * span_ang);
        
        psi = wrap_angle(psi_d + e_psi);
    end
end

function res=sign(a)
    if a >= 0
        res = 1;
    else
        res = -1;

    end
end