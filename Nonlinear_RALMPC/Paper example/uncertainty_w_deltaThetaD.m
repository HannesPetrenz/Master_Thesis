function w_delta_Theta_D = uncertainty_w_deltaThetaD(eta_t,s,z1,z2)
%uncertainty_w_deltaThetaD
%    w_delta_Theta_D = uncertainty_w_deltaThetaD(ETA_T,S,Z1,Z2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    05-Mar-2025 01:05:19

t2 = z1.*6.131811046861033e+1;
t3 = z2.*6.131811046861033e+1;
t4 = z1.*2.129014281432485e+2;
t6 = z2.*2.266033560822537e+1;
t8 = eta_t.*s.*1.252670917307199e-1;
t7 = -t6;
t9 = t3+t4;
t11 = t2+t6;
t14 = z1.*(t3-t4).*(-1.0./2.0e+1);
t12 = t2+t7;
t13 = (t9.*z1)./2.0e+1;
t15 = (t11.*z2)./2.0e+1;
t16 = (t12.*z2)./2.0e+1;
t18 = t13+t15+1.0e-6;
t17 = -t16;
t20 = sqrt(t18);
t19 = t14+t17+1.0e-6;
t21 = eta_t.*t20;
t22 = sqrt(t19);
t24 = t8+t21+1.790989923443473e-5;
t23 = eta_t.*t22;
t25 = t8+t23+1.790989923443473e-5;
w_delta_Theta_D = [t24;t25;t24;t25];
end
