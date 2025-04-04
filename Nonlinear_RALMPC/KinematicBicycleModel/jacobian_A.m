function A = jacobian_A(delta_k,e_psik,e_yk,kappa,v_k)
%jacobian_A
%    A = jacobian_A(DELTA_K,E_PSIK,E_YK,KAPPA,V_K)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    04-Feb-2025 17:51:03

t2 = tan(delta_k);
t3 = e_yk.*kappa;
t4 = t2.^2;
t5 = t2./2.0;
t6 = t3-1.0;
t7 = t4+1.0;
t8 = atan(t5);
t9 = t4./4.0;
t11 = 1.0./t6;
t10 = e_psik+t8;
t12 = t11.^2;
t15 = t9+1.0;
t13 = cos(t10);
t14 = sin(t10);
t16 = 1.0./t15;
t17 = 1.0./sqrt(t15);
A = reshape([1.0,0.0,0.0,0.0,0.0,(kappa.*t12.*t13.*v_k)./1.0e+2,1.0,kappa.^2.*t12.*t13.*v_k.*(-1.0./1.0e+2),0.0,0.0,(t11.*t14.*v_k)./1.0e+2,(t13.*v_k)./1.0e+2,kappa.*t11.*t14.*v_k.*(-1.0./1.0e+2)+1.0,0.0,0.0,t11.*t13.*(-1.0./1.0e+2),t14./1.0e+2,(t2.*t17)./2.5e+1+(kappa.*t11.*t13)./1.0e+2,1.0,0.0,(t7.*t11.*t14.*t16.*v_k)./2.0e+2,(t7.*t13.*t16.*v_k)./2.0e+2,(t7.*t17.*v_k)./2.5e+1-(t4.*t7.*t17.^3.*v_k)./1.0e+2-(kappa.*t7.*t11.*t14.*t16.*v_k)./2.0e+2,0.0,1.0],[5,5]);
end
