% % PT=zeros(10,10,26);
% 

% s=10;
% %Z_10_all=[Z1(:,s); Z2(:,s); Z3(:,s); Z4(:,s); Z5(:,s); Z6(:,s); Z7(:,s); Z8(:,s); Z9(:,s); Z10(:,s); Z11(:,s); Z12(:,s); Z13(:,s); Z14(:,s); Z15(:,s)];
% Z_10_all=[Z1(:,:) Z2(:,:) Z3(:,:) Z4(:,:) Z5(:,:) Z6(:,:) Z7(:,:) Z8(:,:) Z9(:,:) Z10(:,:) Z11(:,:) Z12(:,:) Z13(:,:) Z14(:,:) Z15(:,:)];


% color=['r','g','b','m','c','r','m','g','c','r','b','m','r','g','b','m'];
% color=hsv(16);
% syms x;
% figure;
% j=1;
% clear drift;
% for i=1:2:20
%   y=Z_10_all(i,281:300);
%   [xout b]=hist(y,-1.4:0.2:1.4);
%   summa=sqrt(sum(y.*y));
%   sigma=std(y);
%   avg=mean(y);
%   A=double(int(exp(-(x^2)/(2*sigma^2)),-inf,inf));
%   drift(j,1)=avg;
%   drift(j,2)=2/A+4*i;
%   xx=-1.4:0.01:1.4;
%   c=color(j,:);
%   plot(xx,4*i+2/A*exp((-(xx-avg).^2)/(2*sigma^2)),'Color',c,'LineWidth',1.2,'LineStyle','-.'); hold on;
%   plot(b,xout/summa+4*(i),'Color',c,'LineWidth',1.2); hold on;
%   j=j+1;
% end
% 
% plot([0,0],[0,85],'b','LineWidth',1); hold on;
% plot(drift(:,1),drift(:,2),'--','LineWidth',1);
% title('Ensemble Z15')

% clear PTk;
% j=1;
% for i=150000:1000:160000
%    PTk(:,:,j)=hardy([20 20 1],sdata,data1,i);
%    j=j+1
% end
% 
% X=0:0.003*1000:0.003*10000;
% for i=1:20
%     for j=1:20
%         Y=PTk(i,j,:);
%         Z19(i,j)=trapz(X,Y)/(0.003*10000);
%     end
% end


 
clear gstat;
for k=1:20
    syms x;
    sigma=std(Z19(k,:));
    avg=mean(Z19(k,:));
    A=double(int(exp(-(x^2)/(2*sigma^2)),-inf,inf));
    gstat(k,1)=A;     %save normalization
    gstat(k,2)=avg;   %save average
    gstat(k,3)=sigma;  %save rms
end


color=['r','g','b','m','c','r','m','g','c','r','b','m','r','g','b','m'];
color=hsv(7);
figure; 
j=1;
for k=1:4:20
    x=-5*gstat(k,3):0.001:5*gstat(k,3);
    plot(x,1/gstat(k,1)*exp((-(x-gstat(k,2)).^2)/(2*gstat(k,3)^2)),'Color',color(j,:));
    hold on;
    j=j+1;
end

legend('y=0.1','y=0.3','y=0.5','y=0.7','y=0.9');
title('10000 timesteps integration');
