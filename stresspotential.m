function Sv=stresspotential(Fx,Fy,Fz,xij,yij,zij,lam,vol)
Sv=zeros(3,3);

% Sv(1,1)=-0.5*sum(sum(sum(Fx.*xij))); %xx
% Sv(1,2)=-0.5*sum(sum(sum(Fx.*yij))); %xy
% Sv(1,3)=-0.5*sum(sum(sum(Fx.*zij))); %xz
% 
% Sv(2,1)=-0.5*sum(sum(sum(Fy.*xij))); %yx
% Sv(2,2)=-0.5*sum(sum(sum(Fy.*yij))); %yy
% Sv(2,3)=-0.5*sum(sum(sum(Fy.*zij))); %yz
% 
% Sv(3,1)=-0.5*sum(sum(sum(Fz.*xij))); %zx
% Sv(3,2)=-0.5*sum(sum(sum(Fz.*yij))); %zy
% Sv(3,3)=-0.5*sum(sum(sum(Fz.*zij))); %zz

Sv(1,1)=-0.5*sum(sum(sum(Fx.*xij.*lam))); %xx
Sv(1,2)=-0.5*sum(sum(sum(Fx.*yij.*lam))); %xy
Sv(1,3)=-0.5*sum(sum(sum(Fx.*zij.*lam))); %xz

Sv(2,1)=-0.5*sum(sum(sum(Fy.*xij.*lam))); %yx
Sv(2,2)=-0.5*sum(sum(sum(Fy.*yij.*lam))); %yy
Sv(2,3)=-0.5*sum(sum(sum(Fy.*zij.*lam))); %yz

Sv(3,1)=-0.5*sum(sum(sum(Fz.*xij.*lam))); %zx
Sv(3,2)=-0.5*sum(sum(sum(Fz.*yij.*lam))); %zy
Sv(3,3)=-0.5*sum(sum(sum(Fz.*zij.*lam))); %zz

Sv=Sv/vol;

end