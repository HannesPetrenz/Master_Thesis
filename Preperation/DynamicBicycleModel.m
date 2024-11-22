clear
clc
%% Define symbolic variables
syms x1 x2 x3 x4 x5 x6
syms x1_dot x2_dot x3_dot x4_dot x5_dot x6_dot
syms u1 u2
syms kappa I_z l_f m g mu l_r alpha_f alpha_r fp
%% Systemdefinition
x1_dot=(x5*cos(x3)-x6*sin(x3))/(1-x2*kappa);
x2_dot=x5*sin(x3)+x6*cos(x3);
x3_dot=x4-kappa*(x5*cos(x3)-x6*sin(x3))/(1-x2*kappa);
x4_dot=1/I_z*(-l_f*1/2*m*g*mu*fp+l_r*1/2*m*g*mu*fp);
x5_dot=u1+x4*x6;
x6_dot=1/m*(-1/2*m*g*mu*fp*cos(u2)-1/2*m*g*mu*fp)-x4*x5;

f=[x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot];
u=[u1;u2];
x=[x1;x2;x3;x4;x5;x6];
%% Jacobi
A=simplify(jacobian(f,x))