% function a=localizator(r,dlo, dhi)       
% a=0;
% end
% 
% 
% function b=bond_matrix(b)
% b=0;
% end
% 

%---main---------------
% darg [nx, ny, nz] - number of grid lines per dimension
% sdata [3x2] - box dimensions extracted from dump file
% data [....] trimmed lammps dump output without timestamp headers
% N - timestep you want to Hardy
% 
function output=hardy(darg,sdata,data,N)

% output=['1:i(x) ','2:xmin ','3:xmax ',...
%         '4:j(y) ','5:ymin ','6:ymax ',...
%         '7:k(z) ','8:zmin ','9:zmax ',...
%         '10:rho ','11:momx ','12:momy ','13:momz ',...
%         '14:Q1 ','15:Q2 ','16:Q3 ','17:Q4 ','18:Q5 ',...
%         '19:S1 ','20:S2 ','21:S3 ','22:S4 ','23:S5'];

%wee need to compensate atom coordinates because of removal of walls
% yclo=min(data(:,3));
% ychi=max(data(:,3));
yclo=0.06675+0.00001;
ychi=0.9166-0.00001;
dely=ychi-yclo;


%Recalculate rc from lj units to box normal dimension
rc=1.12246;%potential cutoff distance
rcx=rc/(sdata(1,2)-sdata(1,1));
rcy=rc/(sdata(2,2)-sdata(2,1));
rcz=rc/(sdata(3,2)-sdata(3,1));

%data preprocessing
%pad data with tail and head in x direction
%1***************2%     %2***************1%
%1***************2%     %2***************1%
%1***************2% --> %2***************1%
%1***************2%     %2***************1%
%!***************2%     %2***************1%
%---------------------------
pdata=data(find((data(:,12)==N)),:);
head=pdata(find((1-rcx<=pdata(:,2))&(pdata(:,2)<=1)),:);
tail=pdata(find((0<=pdata(:,2))&(pdata(:,2)<=rcx)),:);
%now fix only X variable in normed and real coordinates by:
%substructing/adding 1 from normed and sdata(1,2)/sdata(1,1) from 
%columns 2 and 8

head(:,2)=head(:,2)-1;
head(:,8)=head(:,8)-(sdata(1,2)-sdata(1,1));
tail(:,2)=tail(:,2)+1;
tail(:,8)=tail(:,8)+(sdata(1,2)-sdata(1,1));

pdata=[pdata;head;tail];
%---------------------------


[X Y]=meshgrid(1:darg(1),1:darg(2));
zmesh=meshgrid(1:darg(1),1:darg(2));


%figure;

%output=zeros(darg(1)*darg(2)*darg(3),23);
ii=1;%variable to loop over output
for i=1:darg(1)
    for j=1:darg(2)
        for k=1:darg(3)
            %scaled box dimenstions
            sxlo=(j-1)/darg(1);
            sylo=(i-1)*dely/darg(2)+yclo;
            szlo=(k-1)/darg(3);
            sxhi=(j)/darg(1);
            syhi=(i)*dely/darg(2)+yclo;
            szhi=(k)/darg(3);
            
            dx=(sdata(1,2)-sdata(1,1))*(sxhi-sxlo);
            dy=(sdata(2,2)-sdata(2,1))*(syhi-sylo);
            dz=(sdata(3,2)-sdata(3,1))*(szhi-szlo);
            
            %true box dimensions
            xlo=sdata(1,1)+dx*(j-1);
            xhi=sdata(1,1)+dx*j;
            ylo=sdata(2,1)+dy*(i-1)+(sdata(2,2)-sdata(2,1))*yclo;
            yhi=sdata(2,1)+dy*i+(sdata(2,2)-sdata(2,1))*yclo;
            zlo=sdata(3,1)+dz*(k-1);
            zhi=sdata(3,1)+dz*k;            
            
            dsx=(sxhi-sxlo);
            dsy=(syhi-sylo);
            dsz=(szhi-szlo);
            

            SR=sqrt(dsx^2+dsy^2); %Scaled Radius
            R=sqrt(dx^2+dy^2); %True Radius

            
%            vol=dz*(pi*R^2);
             
            
            vol=(sdata(1,2)-sdata(1,1))*(sxhi-sxlo)*...
                (sdata(2,2)-sdata(2,1))*(syhi-sylo)*...
                (sdata(3,2)-sdata(3,1))*(szhi-szlo);
           
%            %extract atoms at within local cylindrical region
%             psiatoms=data(find((data(:,12)==N)&...
%                       (sqrt((data(:,2)-sxlo).^2+(data(:,3)-sylo).^2)<=SR)),:);


%            %extract atoms within local cubical region
            psiatoms=pdata(find((sxlo<=pdata(:,2))&(pdata(:,2)<=sxhi)&...
                                (sylo<=pdata(:,3))&(pdata(:,3)<=syhi)&...
                                (szlo<=pdata(:,4))&(pdata(:,4)<=szhi)),:);
       
            
%           %extrat atoms within potential cutoff outside of cubical region
            outsiders=pdata(find((sxlo-rcx<=pdata(:,2))&(pdata(:,2)<=sxhi+rcx)&...
                            (sylo-rcy<=pdata(:,3))&(pdata(:,3)<=syhi+rcy)&...
                            (szlo-rcz<=pdata(:,4))&(pdata(:,4)<=szhi+rcz)),:);
            outsiders(ismember(outsiders(:,1),psiatoms(:,1)),:)=[];                        
            
%             outsiders=data(find((data(:,12)==N)&...
%                 ((sxhi<=data(:,2))&(data(:,2)<=sxhi+rcx))|...
%                 ((syhi<=data(:,3))&(data(:,3)<=syhi+rcy))|...
%                 ((szhi<=data(:,4))&(data(:,4)<=szhi+rcz))|...
%                 ((sxlo-rcx<=data(:,2))&(data(:,2)<=sxlo))|...
%                 ((sylo-rcy<=data(:,3))&(data(:,3)<=sylo))|...
%                 ((szlo-rcz<=data(:,4))&(data(:,4)<=szlo))),:);
                


                
           %find potentials forces and distance matrix
           [phi Fx Fy Fz xij yij zij lam]=neighlist(psiatoms,outsiders,[xlo,xhi;ylo,yhi;zlo,zhi]);            
            
            %calculate total values of [mass vx vy vz] for 2d box
            totals=sum(psiatoms(:,[11 5 6 7]));

           
%              %calculate total values of [mass vx vy vz] for box
%              totals=sum(data(find((data(:,12)==N)&...
%                     (sxlo<=data(:,2))&(data(:,2)<=sxhi)&...
%                     (sylo<=data(:,3))&(data(:,3)<=syhi)&...
%                     (szlo<=data(:,4))&(data(:,4)<=szhi)),...
%                                                   [11 5 6 7]));
             
             %calculate number of particles
             NN=length(psiatoms(:,1));      
             
             %separate calculated values
             totmass=totals(1);
             avevelx=totals(2)/NN;
             avevely=totals(3)/NN;
             avevelz=totals(4)/NN;

           %find density
           rho=(totmass/vol);
           
           %find momentum
           momx=totmass*avevelx/(vol);
           momy=totmass*avevely/(vol);
           momz=totmass*avevelz/(vol);
           
%            Qk=heatkinetic(psiatoms,[avevelx avevely avevelz],phi,vol);
%            
%            Qv=heatpotential(psiatoms,[avevelx avevely avevelz],Fx,Fy,Fz,xij,yij,zij,lam,vol);
           
           Sk=stresskinetic(psiatoms,[avevelx avevely avevelz],vol);
           
           Sv=stresspotential(Fx,Fy,Fz,xij,yij,zij,lam,vol);
           
           P0=-trace(Sk+Sv)/3;
           
           %find heat tensor
           Q=zeros(1,5);
           
           %find stress tensor
           S=zeros(1,5);
             
           SS=(-Sv-Sk);
           zmesh(i,j)=SS(1,2);
           
           if isnan(zmesh(i,j))
               yclo=yclo-yclo+yclo; %catch NaN
           end
           
           %output(i,j)=trace(SS)/3;
           
%            output(ii,:)=[i, 0, 0, j, 0, 0, k, 0, 0,...
%                         rho, momx, momy, momz, Q, S];
                     
                    
%----------------------------
% %            plot(data(find((data(:,12)==N)&...
% %                     (sxlo<=data(:,2))&(data(:,2)<=sxhi)&...
% %                     (sylo<=data(:,3))&(data(:,3)<=syhi)&...
% %                     (szlo<=data(:,4))&(data(:,4)<=szhi)),...
% %                                                   2),...
% %                     data(find((data(:,12)==N)&...
% %                     (sxlo<=data(:,2))&(data(:,2)<=sxhi)&...
% %                     (sylo<=data(:,3))&(data(:,3)<=syhi)&...
% %                     (szlo<=data(:,4))&(data(:,4)<=szhi)),...
% %                                                   3),...
% %                     'LineStyle','.','Color',[rand() rand() rand()]);
% %            hold on;
%----------------------------
                    
                    
           ii=ii+1;
        end
    end
end

% figure;
% surf(X,Y,zmesh);

output=zmesh;

end
