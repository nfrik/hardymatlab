function Sk=stresskinetic(psiatoms,vel,vol)
Sk=zeros(3,3);
v=[psiatoms(:,5) psiatoms(:,6) psiatoms(:,7)];
m=psiatoms(:,11);
vminusu=[v(:,1)-vel(1) v(:,2)-vel(2) v(:,3)-vel(3)];
for i=1:3
    for j=1:3
        Sk(i,j)=-sum(m.*vminusu(:,i).*vminusu(:,j));
    end
end
Sk=Sk/vol;

end