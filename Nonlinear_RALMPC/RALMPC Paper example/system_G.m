function G = system_G(x_1k,x_2k)
%system_G
%    G = system_G(X_1K,X_2K)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    06-Mar-2025 20:56:08

G = reshape([x_2k.*(-1.0./2.0e+1),0.0,0.0,x_1k./2.0e+1],[2,2]);
end
