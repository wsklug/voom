clear;

dt = 1.0e-4;
kT=4.1;
gamma=0.6;

load txy.dat;

sz=size(txy);
npts=(sz(2)-1)/2;
n=sz(1);

t=dt*txy(:,1);
x=txy(:,2:2:sz(2));
y=txy(:,3:2:sz(2));
r2=x.^2+y.^2;
r=sqrt(r2);



% a=0;
% b=1;
% dz = a + (b-a).*rand(n,1);
% z = 0;
% for i=1:n,
%     if(dz(i)<0.5) z = [z, z(i)-1];
%     else z = [z, z(i)+1];
%     end
% end
Xavg=zeros(n,1);
Yavg=zeros(n,1);
Ravg=zeros(n,1);
dX2=zeros(n,1);
dY2=zeros(n,1);
dR2=zeros(n,1);

for i=1:n,
    Xavg(i)=sum(x(i,:))/npts;
    Yavg(i)=sum(y(i,:))/npts;
    Ravg(i)=sum(r(i,:))/npts;
    dX2(i)=sum( x(i,:).^2 )/npts;
    dY2(i)=sum( y(i,:).^2 )/npts;
    dR2(i)=sum( r(i,:).^2 )/npts;
end

figure(1);
subplot(311);
loglog(t,dX2,'.-');
subplot(312);
loglog(t,dY2,'.-');
subplot(313);
loglog(t,dR2,'.-');
    
figure(2);
plot(x(:,1),y(:,1),'-');
hold on;
plot(x(1,1),y(1,1),'g*');
plot(x(n,1),y(n,1),'r*');
axis equal;
hold off;

figure(3);
subplot(311);
plot(t,x(:,1),'.-');
subplot(312);
plot(t,y(:,1),'.-');
subplot(313);
plot(t,r(:,1),'.-');

figure(4);
subplot(311);
plot(t,dX2,'.-',t,2*kT/gamma*t,'-.');
subplot(312);
plot(t,dY2,'.-',t,2*kT/gamma*t,'-.');
subplot(313);
plot(t,dR2,'.-',t,4*kT/gamma*t,'-.');

figure(5);
subplot(311);
plot(t,Xavg,'.-');
subplot(312);
plot(t,Yavg,'.-');
subplot(313);
plot(t,Ravg,'.-');
