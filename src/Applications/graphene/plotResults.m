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
plot(z,s0max,'b.-',z,s1max,'b--');
hold on;

figure(3);
plot(-load,s0max,'b.-',-load,s1max,'b--');
hold on;

clearvars zfs step z load restraint s0min s0max s1min s1max;

zfs=load('afmR150/afmR150.zfs');
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
plot(z,s0max,'r.-',z,s1max,'r--');

figure(3);
plot(-load,s0max,'r.-',-load,s1max,'r--');

clearvars zfs step z load restraint s0min s0max s1min s1max;

zfs=load('afmR200/afmR200.zfs');
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
legend('R=172nm','R=150nm','R=200nm');

figure(2);
plot(z,s0max,'g.-',z,s1max,'g--');
hold off;

axis([0,zmax,0,3.5e4]);
legend('zz max R=172nm', 'ac max R=172nm', 'zz max R=150nm', 'ac max R=150nm', 'zz max R=200nm', 'ac max R=200nm');
xlabel('Indentation Depth [nm]');
ylabel('Stress [pN/nm]');

figure(3);
plot(-load,s0max,'g.-',-load,s1max,'g--');
hold off;

axis([0,16e6,0,3.5e4]);
legend('zz max R=172nm', 'ac max R=172nm', 'zz max R=150nm', 'ac max R=150nm', 'zz max R=200nm', 'ac max R=200nm');
xlabel('Load [pN]');
ylabel('Stress [pN/nm]');

