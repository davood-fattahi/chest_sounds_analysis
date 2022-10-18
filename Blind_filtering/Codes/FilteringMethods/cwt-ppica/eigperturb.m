function [V,D]=eigperturb(V0,D0,dA,dB)
% 
% update the generalized eigen decomposition with perturbation in 
% the matrices. The following function is based on perturbation theory 
% formulation.
%
% Inputs
% V0: the initial eigenvectors
% D0: the initial eigenvalues
% dA: perturbation in matrix A
% dB: perturbation in matrix B
% 
% Outputs
% V: the updated eigenvectors
% D: the updated eigenvalues
% 
% 
% Davood Fattahi, fattahi.d@gmail.com, 14/05/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isdiag(D0)
    D0=diag(D0);
end
I=1:size(D0,1);
for i=I
    D(i)=D0(i)+V0(:,i)'*(dA-D0(i)*dB)*V0(:,i);
    S=0;
    for j=I(I~=i)
        S=S + (1/(D0(i)-D0(j)))*(V0(:,j)'*(dA-D0(i)*dB)*V0(:,i)*V0(:,j));
    end
    V(:,i)=V0(:,i)*(1-.5*V0(:,i)'*dB*V0(:,i))+S;
end

D=diag(D);