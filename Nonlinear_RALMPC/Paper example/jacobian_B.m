function B = jacobian_B(x_1k,x_2k)
%jacobian_B
%    B = jacobian_B(X_1K,X_2K)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    06-Mar-2025 20:56:08

B = [x_1k./4.0e+1+1.0./4.0e+1;x_2k.*(-1.0./1.0e+1)+1.0./4.0e+1];
end
