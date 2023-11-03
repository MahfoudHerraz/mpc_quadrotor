function Hess_J = HessianJ(x,lambda,Pm,Psi_map,Cd_map)
Hess_J = 2*Pm'*Psi_map'*(Cd_map'*Cd_map)*Psi_map*Pm ...
    + 2*(Pm'*Pm);
end