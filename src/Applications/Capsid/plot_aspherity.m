%
%
%
sub4=load('hk97/relax/sub4.aspherity.dat');
g=sub4(:,1);
a=sub4(:,2);
%
%
%
sub4w=load('t4/relax/sub4w.aspherity.dat');
gT4w=sub4w(:,1);
aT4w=sub4w(:,2);
%
%
%
sub4w=load('hk97/relax/sub4w.aspherity.dat');
gT7w=sub4w(:,1);
aT7w=sub4w(:,2);
%
%
%
sub3w=load('t16/relax/sub3w.aspherity.dat');
gT16w3=sub3w(:,1);
aT16w3=sub3w(:,2);
%
%
%
LW = 1.5;
MS = 6;
semilogx(g,a,'ko-',  'LineWidth', LW, 'MarkerSize', MS, 'MarkerFaceColor', 'w');
hold on;

semilogx(gT16w3,aT16w3,'ro-',  'LineWidth', LW, 'MarkerSize', MS, 'MarkerFaceColor', 'r');
semilogx(gT7w,aT7w,'b^-',  'LineWidth', LW, 'MarkerSize', MS, 'MarkerFaceColor', 'b');
semilogx(gT4w,aT4w,'gs-',  'LineWidth', LW, 'MarkerSize', MS, 'MarkerFaceColor', 'g');
legend('Intact','T=16', 'T=7', 'T=4');
% ylabel('< \Delta R^2 >/ < R^2 >');
% xlabel('\gamma = YR^2/\kappa');
ylabel('DR');
xlabel('g');
axis([1 10^8 0 0.002]);
hold off;