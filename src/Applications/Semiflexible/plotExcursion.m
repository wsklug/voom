clear;

dt = 1.0e-4;
kT=4.1;
gamma=0.6;

load txy.dat;

t=dt*txy(:,1);
x=txy(:,2);
y=txy(:,3);
r2=x.^2+y.^2;
r=sqrt(r2);

n=size(t);
n=n(1);

% a=0;
% b=1;
% dz = a + (b-a).*rand(n,1);
% z = 0;
% for i=1:n,
%     if(dz(i)<0.5) z = [z, z(i)-1];
%     else z = [z, z(i)+1];
%     end
% end

nplot = 1000;
n = n/nplot;
for i=1:nplot,
    T(i)=t(n*i);
    Xavg(i)=sum(x(1:n*i))/(n*i);
    Yavg(i)=sum(y(1:n*i))/(n*i);
%     Zavg(i)=sum(z(1:n*i))/(n*i);
    Ravg(i)=sum(r(1:n*i))/(n*i);
%    dX2(i)=sum( (x(1:n*i)-Xavg(i)).^2 );
%    dY2(i)=sum( (y(1:n*i)-Yavg(i)).^2 );
%    dR2(i)=sum( (r(1:n*i)-Ravg(i)).^2 );
    dX2(i)=sum( (x(1:n*i)).^2 )/(n*i);
    dY2(i)=sum( (y(1:n*i)).^2 )/(n*i);
%     dZ2(i)=sum( (z(1:n*i)).^2 )/(n*i);
    dR2(i)=sum( (r(1:n*i)).^2 )/(n*i);
end

figure(1);
subplot(411);
loglog(T,dX2,'.-');
subplot(412);
loglog(T,dY2,'.-');
subplot(413);
loglog(T,dR2,'.-');
% subplot(414);
% loglog(T,dZ2,'.-');
    
figure(2);
plot(x,y,'-');
hold on;
plot(x(1),y(1),'g*');
plot(x(n*nplot),y(n*nplot),'r*');
axis equal;
hold off;

figure(3);
subplot(311);
plot(t,x,'.-');
subplot(312);
plot(t,y,'.-');
subplot(313);
plot(t,r,'.-');

figure(4);
subplot(311);
plot(T,dX2,'.-',T,2*kT/gamma*T,'-.');
subplot(312);
plot(T,dY2,'.-',T,2*kT/gamma*T,'-.');
subplot(313);
plot(T,dR2,'.-',T,4*kT/gamma*T,'-.');

