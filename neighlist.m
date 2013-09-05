% neighlist generates matrix of potentials and can be used as a bond matrix

function [u,fx,fy,fz,rijx,rijy,rijz,lam]=neighlist(psiatoms,outsiders,boxdim)
rowspsi=length(psiatoms(:,1));
rowsout=length(outsiders(:,1));
rows=rowspsi+rowsout;
u=zeros(rows,rows);
fx=zeros(rows,rows);
fy=zeros(rows,rows);
fz=zeros(rows,rows);
rijx=zeros(rows,rows);
rijy=zeros(rows,rows);
rijz=zeros(rows,rows);
lam=zeros(rows,rows);

%figure; plot(outsiders(:,8),outsiders(:,9),'.','Color',[rand() rand() rand()]); hold on; plot(psiatoms(:,8),psiatoms(:,9),'*','Color',[rand() rand() rand()]);

for i=1:rowspsi
    for j=i+1:rowspsi
        if i~=j            
            %A=[psiatoms(i,9) psiatoms(i,10) psiatoms(i,11)];
            %B=[psiatoms(j,9) psiatoms(j,10) psiatoms(j,11)];
            A=psiatoms(i,:);
            B=psiatoms(j,:);
            u(i,j)=potential(A(8:10),B(8:10));            
            [fx(i,j) fy(i,j) fz(i,j)]=force(A(8:10),B(8:10));
            lam(i,j)=1;
            r=A(8:10)-B(8:10);
            rijx(i,j)=r(1);
            rijy(i,j)=r(2);
            rijz(i,j)=r(3);
        end
    end
end


for i=1:rowspsi
    for j=1:rowsout
            %A=[psiatoms(i,9) psiatoms(i,10) psiatoms(i,11)];
            %B=[psiatoms(j,9) psiatoms(j,10) psiatoms(j,11)];
            l=i;
            m=j+rowspsi;
            A=psiatoms(i,:);
            B=outsiders(j,:);
            u(l,m)=potential(A(8:10),B(8:10));            
            [fx(l,m) fy(l,m) fz(l,m)]=force(A(8:10),B(8:10));
            if (fx(l,m)~=0)||(fy(l,m)~=0)||(fz(l,m)>0)
            [lam(l,m) mx my mz]=lambda(A(8:10),B(8:10),boxdim);%need to pass scaled coords
            %lam(l,m)=0;
            %plot([A(8) B(8)],[A(9) B(9)],'-','Color','b');
            %plot([A(8) mx],[A(9) my],'-','Color','r');
            end
            r=A(8:10)-B(8:10);
            rijx(l,m)=r(1);
            rijy(l,m)=r(2);
            rijz(l,m)=r(3);
    end
end

u=u+u'; %potential is symmetric
fx=fx-fx'; %force is antysymmetric
fy=fy-fy'; %force is antysymmetric
fz=fz-fz'; %force is antysymmetric
rijx=rijx-rijx'; %distance is antysymmetric
rijy=rijy-rijy'; %distance is antysymmetric
rijz=rijz-rijz'; %distance is antysymmetric
lam=lam+lam';
end