%--------------------------
%A/B - position of the particle inside/outside the box
%boxdim - dimension of the box
%lam=|A-M|/|A-B| where M is the point laying on the border of the box
function [lam mx my mz]=lambda(A,B,boxdim)
xlo=boxdim(1,1);
xhi=boxdim(1,2);
ylo=boxdim(2,1);
yhi=boxdim(2,2);
% zlo=boxdim(3,1);
% zhi=boxdim(3,2);

k=(B(2)-A(2))/(B(1)-A(1));
b=B(2)-k*B(1);

ym1=k*xlo+b;
xm1=xlo;
ym2=ylo;
xm2=(ylo-b)/k;
ym3=yhi;
xm3=(yhi-b)/k;
ym4=k*xhi+b;
xm4=xhi;

%plot([boxdim(1,1) boxdim(1,2) boxdim(1,2) boxdim(1,1) boxdim(1,1)],[boxdim(2,1) boxdim(2,1) boxdim(2,2) boxdim(2,2) boxdim(2,1)],'-','Color','b'); hold on;

M=[xm1 ym1 0;xm2 ym2 0;xm3 ym3 0;xm4 ym4 0];

for i=1:4
    Mm(i)=sum(abs(B-M(i,:))+abs(A-M(i,:)));
end


% if (ym1<=yhi)&&(ym1>=ylo)&&(xm1<=xhi)&&(xm1>=xlo)
%    M=[xm1, ym1, 0];
%    lamb(1)=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
% elseif (ym2<=yhi)&&(ym2>=ylo)&&(xm2<=xhi)&&(xm2>=xlo)
%    M=[xm2, ym2, 0];
%    lamb=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
%    lamb(2)=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
% elseif (ym3<=yhi)&&(ym3>=ylo)&&(xm3<=xhi)&&(xm3>=xlo)
%    M=[xm3, ym3, 0];
%    lamb=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
%    lamb(3)=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
% elseif (ym4<=yhi)&&(ym4>=ylo)&&(xm4<=xhi)&&(xm4>=xlo) 
%    M=[xm4, ym4, 0];
%    lamb=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
%    lamb(4)=sqrt(sum((A-M).*(A-M)))/sqrt(sum((A-B).*(A-B)));
% end

[lam indx]=min(Mm);
lam=sqrt(sum((A-M(indx,:)).*(A-M(indx,:))))/sqrt(sum((A-B).*(A-B)));
% if ~isreal(lam)
%     lam=lam'';
% end
mx=M(indx,1);
my=M(indx,2);
mz=M(indx,3);
%plot(mx,my,'.','Color','g');
end