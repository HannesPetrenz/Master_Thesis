function B = jacobian_B(theta)
%jacobian_B
%    B = jacobian_B(THETA)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    04-Feb-2025 17:51:03

B = reshape([0.0,0.0,0.0,0.0,1.0./1.0e+2,0.0,0.0,0.0,theta./1.0e+2,0.0],[5,2]);
end
