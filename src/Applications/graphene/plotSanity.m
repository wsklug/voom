clearvars zfs step z load restraint s0min s0max s1min s1max;

zfs=load('hex30Q5thOrder/hex30Q5thOrder.zfs');
step=zfs(:,1);
z=zfs(:,2);
load=zfs(:,3);
restraint=zfs(:,4);
s0min=zfs(:,5);
s0max=zfs(:,6);
s1min=zfs(:,7);
s1max=zfs(:,8);

figure(1);
plot(z,-load,'b.-');
hold on;

figure(2);
plot(z,s0max,'b.-'); %,z,s1max,'b--');
hold on;

figure(3);
plot(-load,s0max,'b.-');%,-load,s1max,'b--');
hold on;

clearvars zfs step z load restraint s0min s0max s1min s1max;

%zfs=load('r175p7e-4/hex30Q.zfs');
zfs=load('r105/hex30Q.zfs');
step=zfs(:,1);
z=zfs(:,2);
load=zfs(:,3);
restraint=zfs(:,4);
s0min=zfs(:,5);
s0max=zfs(:,6);
s1min=zfs(:,7);
s1max=zfs(:,8);

figure(1);
plot(z,-load,'r.-');
hold on;

figure(2);
plot(z,s0max,'r.-'); %,z,s1max,'r--');

figure(3);
plot(-load,s0max,'r.-'); %,-load,s1max,'r--');

clearvars zfs step z load restraint s0min s0max s1min s1max;

%zfs=load('r175p3e-3/hex30Q.zfs');
zfs=load('r27-5/hex30Q.zfs');
step=zfs(:,1);
z=zfs(:,2);
load=zfs(:,3);
restraint=zfs(:,4);
s0min=zfs(:,5);
s0max=zfs(:,6);
s1min=zfs(:,7);
s1max=zfs(:,8);

figure(1);
plot(z,-load,'g.-');
hold off;

zmax = 220;

axis([0,zmax,0,16e6]);
xlabel('Indentation Depth [nm]');
ylabel('Load [pN]');
%legend('\tau=0 R=172nm','\tau=1.5 N/m R=130nm','\tau=1.5 N/m R=200nm');
legend('\tau=0 R=172nm','\tau=1.5 N/m R=105nm','\tau=0.34 N/m R=27.5nm');

figure(2);
plot(z,s0max,'g.-'); %,z,s1max,'g--');
hold off;

axis([0,zmax,0,3.5e4]);
%legend('zz max \tau=0', 'ac max \tau=0', 'zz max \tau=1.5 N/m R=130nm', 'ac max \tau=1.5 N/m R=130nm', 'zz max \tau=1.5 N/m R=200nm', 'ac max \tau=1.5 N/m R=200nm');
%legend('zz max \tau=0 R=172nm',  'zz max \tau=1.5 N/m R=130nm', 'zz max \tau=1.5 N/m R=200nm');
legend('\tau=0 R=172nm','\tau=1.5 N/m R=105nm','\tau=0.34 N/m R=27.5nm');
xlabel('Indentation Depth [nm]');
ylabel('Stress [pN/nm]');

figure(3);
plot(-load,s0max,'g.-'); %,-load,s1max,'g--');
hold off;

axis([0,16e6,0,3.5e4]);
%legend('zz max \tau=0', 'ac max \tau=0', 'zz max \tau=0.4 N/m', 'ac max \tau=0.4 N/m', 'zz max \tau=1.5 N/m', 'ac max \tau=1.5 N/m');
%legend('zz max \tau=0 R=172nm',  'zz max \tau=1.5 N/m R=130nm', 'zz max \tau=1.5 N/m R=200nm');
legend('\tau=0 R=172nm','\tau=1.5 N/m R=105nm','\tau=0.34 N/m R=27.5nm');
xlabel('Load [pN]');
ylabel('Stress [pN/nm]');

