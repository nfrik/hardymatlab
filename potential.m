function u=potential(A,B)
sigma=1.0;
epsilon=1.0;
rc=1.12246;
r=sqrt(sum(((A-B).^2)'));
% u=4*epsilon*((sigma./r).^12-(sigma./r).^6);
% zeros=r<rc;
% u=u.*zeros;
u=0;
if r<rc
   u=4*epsilon*((sigma/r)^12-(sigma/r)^6);
end
end
