function [J,grad_J,Hess_J] = cost_function(Phi_map,Psi_map,Cd_map,p,Yt,X,Pm)
% v = (P*reshape(p,n,4))';
% v = reshape(Pm*p,4,N);
v = Pm*p;
Y = Cd_map*Phi_map*X + Cd_map*Psi_map*v;
J = norm(Y-Yt(:))^2 + norm(v)^2;

if nargout>1
    grad_J = 2*Pm'*Psi_map'*Cd_map'*(Y-Yt(:)) ...
        + 2*Pm'*v;
end

if nargout>2
    Hess_J = 2*Pm'*Psi_map'*(Cd_map'*Cd_map)*Psi_map*Pm ...
        + 2*(Pm'*Pm);
end
% Y = Cd*X;
% Xold = X;
% J = norm(Yt(:,1)-Y)+0*norm(v(:,1));
% for i=2:N
%     Xnew = Ad*Xold+Bd*v(:,i-1);
%     Y = Cd*Xnew;
%     J = J+norm(Yt(:,i)-Y)+0*norm(v(:,i));
%     Xold = Xnew;
% end
end