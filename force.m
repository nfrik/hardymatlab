function [Fx Fy Fz]=force(A,B)
sigma=1.0;
epsilon=1.0;
rc=1.12246;
r=sqrt(sum(((A-B).^2)'));
rr=(A-B)/r;
% u=4*epsilon*((sigma./r).^12-(sigma./r).^6);
% zeros=r<rc;
% u=u.*zeros;
f=[0 0 0];
if r<rc
   %u=4*epsilon*((sigma/r)^12-(sigma/r)^6);
   f=24*epsilon*(2*(sigma^12)/(r^13)-(sigma^6)/(r^7))*rr;
end

Fx=f(1);
Fy=f(2);
Fz=f(3);

end
