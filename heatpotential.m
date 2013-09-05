% psiatoms - atoms in the area
% vel - averaged cell velocity [vx vy vz]
% F - potentials matrix
% rij - interatom distance matrix
% vol - volume
% Qk  - vector [Qk_x Qk_y Qk_z]

function Qv=heatpotential(psiatoms,vel,Fx,Fy,Fz,xij,yij,zij,lam,vol)
Qv_x=0;
Qv_y=0;
Qv_z=0;
Qv=0;
r=[psiatoms(:,8) psiatoms(:,9) psiatoms(:,10)];
v=[psiatoms(:,5) psiatoms(:,6) psiatoms(:,7)];
vminusu=[v(:,1)-vel(1) v(:,2)-vel(2) v(:,3)-vel(3)];

Qv_x=0.5*sum((Fx'.*xij.*lam)*vminusu(:,1)+(Fy'.*xij.*lam)*vminusu(:,2)+(Fz'.*xij.*lam)*vminusu(:,2));
Qv_y=0.5*sum((Fx'.*yij.*lam)*vminusu(:,1)+(Fy'.*yij.*lam)*vminusu(:,2)+(Fz'.*yij.*lam)*vminusu(:,2));
Qv_z=0.5*sum((Fx'.*zij.*lam)*vminusu(:,1)+(Fy'.*zij.*lam)*vminusu(:,2)+(Fz'.*zij.*lam)*vminusu(:,2));

% Qv_x=0.5*sum(Fx*vminusu(:,1))/vol;
% Qv_y=0.5*sum(Fy*vminusu(:,2))/vol;
% Qv_z=0.5*sum(Fz*vminusu(:,3))/vol;

% Qv=0.5*F.*rij;
% Qv_x=sum(Qv*vminusu(:,1)');
% Qv_y=sum(Qv*vminusu(:,2)');
% Qv_z=sum(Qv*vminusu(:,3)');
Qv=[Qv_x Qv_y Qv_z]/vol;
end