% psiatoms - atoms in the area
% vel - average cell velocity [vx vy vz]
% phi - potentials matrix
% vol - volume
% Qk  - vector [Qk_x Qk_y Qk_z]

function Qk=heatkinetic(psiatoms,vel,phi,vol)
Qk_x=0;
Qk_y=0;
Qk_z=0;
Qk=0;
r=[psiatoms(:,8) psiatoms(:,9) psiatoms(:,10)];
v=[psiatoms(:,5) psiatoms(:,6) psiatoms(:,7)];
m=psiatoms(:,11);
vminusu=[v(:,1)-vel(1) v(:,2)-vel(2) v(:,3)-vel(3)];
Qk=(0.5*m.*sum(vminusu.^2')'+0.5*sum(phi')');
Qk_x=Qk.*vminusu(:,1);
Qk_y=Qk.*vminusu(:,2);
Qk_z=Qk.*vminusu(:,3);
Qk=[sum(Qk_x) sum(Qk_y) sum(Qk_z)]/vol;
end