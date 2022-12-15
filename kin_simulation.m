 %% 3dof Robotic Manipulator
clear all;
close all;
l(1) = 5.0;
l(2) = 5.0; %% in cm
l(3) = 5.0;
l(4) = 5.0; 

%% sampling period
dt = 0.001; %dt=1 msec
%%  workspace
Tf=10.0; % movement period
t=0:dt:Tf;
%xdA,xdB,ydA,ydB,zdA,zdB: starting and end point, end-effector
xdA = -20.0;
xdB = 20.0;
ydA = -10.00;
ydB = 10.00;
zdA = 30.0;
zdB = 30.0;

% position management
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...');
disp(' ');
xd(1) = xdA;
yd(1) = ydA;
zd(1) = zdA;
lambda_x = (xdB-xdA)/Tf;
lambda_y = (ydB-ydA)/Tf;
lambda_z = (zdB-zdA)/Tf;
kmax=Tf/dt + 1;
for k=2:kmax;
xd(k) = xd(k-1) + lambda_x*dt; 
yd(k) = yd(k-1) + lambda_y*dt; 
zd(k) = zd(k-1) + lambda_z*dt; 
end 

%% ******Kinematic Simulation ******
disp('Kinematic Simulation ...'); %%
disp(' '); %%


Pez=zdB;
% smooth the motion through polynomial interpolation functions

dPeyT0=0;%dPey(0) /dt =0
dPeyTf=0; %dPey(10)/dt=0

syms a0 a1 a2 a3
a0=ydA; 
a1=dPeyT0;
a2=(3/Tf^2)*(ydB-ydA)-(2/Tf)*a1-(1/Tf)*dPeyTf;
a3=(-2/Tf^3)*(ydB-ydA)+(1/10^2)*(dPeyT0+dPeyTf);

dPexT0=0;
dPexTf=0;

syms bo b1 b2 b3
b0=xdA;
b1=dPexT0;
b2=(3/Tf^2)*(xdB-xdA)-(2/Tf)*b1-(1/Tf)*dPexTf;
b3=(-2/Tf^3)*(xdB-xdA)+(1/10^2)*(dPexT0+dPexTf);

Pex=b0+b1*t(:)+b2*t(:).^2+b3*t(:).^3;
Pey=a0+a1*t(:)+a2*t(:).^2+a3*t(:).^3;

%% ***** Inverse Kinematic Simulation *****
%% joints
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., kmax,}
qd(1:kmax,1)=asin(-l(2)./(sqrt(Pex(:).^2+(l(1)-Pez)^2)))-asin((l(1)-Pez)./(sqrt(Pex(:).^2+(l(1)-Pez)^2)));
s1=sin(qd(1,1));  % s1 c1 constants because q1 is constant (independet of q2 ,q3
c1=cos(qd(1,1)); 
qd(:,3)=acos( (Pey(:).^2+( (Pez-c1*l(2)-l(1))/s1)^2-l(4)^2+l(3)^2 )./ (2*l(4)*l(3)));
s3=sin(qd(:,3));
c3=cos(qd(:,3));
qd(:,2)=asin( (Pey(:)./(sqrt( (l(4).*s3(:)).^2+( l(3)+l(4).*c3(:)).^2))))- ...
       asin( (l(4).*s3(:)./(sqrt( (l(4).*s3(:)).^2+( l(3)+l(4).*c3(:)).^2))));
s2=sin(qd(:,2));
c2=cos(qd(:,2));
%% ***** Kinematic Analysis*****
%%xdXY=cartesian positions of Ax->y(q) matrix
xd00(1:kmax)= 0;
yd00(1:kmax)=0;
zd00(1:kmax)=0;
xd01(1:kmax)= 0;
yd01(1:kmax) =0;
zd01(1:kmax)=l(1);
xd02(1:kmax) =-l(2)*s1;
yd02(1:kmax) =0;
zd02(1:kmax)=l(2)*c1+l(1);
xd03 = c1.*c2(:).*l(3)-l(2)*s1 ;
yd03 = l(3).*s2(:)
zd03= l(3)*s1.*c2(:)+l(2)*c1+l(1);
c23=cos(qd(:,2)+qd(:,3));
s23=sin(qd(:,2)+qd(:,3));
xd0E=c1.*c23(:).*l(4)+c1.*c2(:).*l(3)-s1*l(2)
yd0E=s23(:).*l(4)+s2(:).*l(3);
zd0E=s1.*c23(:).*l(4)+s1.*c2(:).*l(3)+c1*l(2)+l(1); 
%% *** PLOTS *** %%**

fig1 = figure;
subplot(2,2,1);
plot(t,xd);
ylabel('xd (cm)');
xlabel('time t (sec)');

subplot(2,2,2);
plot(t,yd);
ylabel('yd (cm)');
xlabel('time t (sec)');

subplot(2,2,3);
plot(t,zd);
ylabel('zd (cm)');
xlabel('time t (sec)');

fig2=figure;
Vx=diff(Pex);
Vx=[0;Vx];
Vy=diff(Pey);
Vy=[0;Vy];
Vz=diff(Pez);
Vz=[0;Vz];
subplot(2,2,1);
plot(t,Vx);
ylabel('Vx (cm/sec)');
xlabel('time t (sec)');
subplot(2,2,2);
plot(t,Vy);
ylabel('Vy (cm/sec)');
xlabel('time t (sec)');
subplot(2,2,3);
plot(t,Vz);
ylabel('Vz (cm/sec)');
xlabel('time t (sec)'); 

fig3 = figure;
subplot(2,2,1);
plot(t,(qd(:,1)));
ylabel('qd1 (rad)');
xlabel('time t (sec)');

subplot(2,2,2);
plot(t,qd(:,2));
ylabel('qd2 (rad)');
xlabel('time t (sec)');

subplot(2,2,3);
plot(t,qd(:,3));
ylabel('qd3 (rad)');
xlabel('time t (sec)');


fig4 = figure;
dif_qd1=diff(qd(:,1));
dif_qd2=diff(qd(:,2));
dif_qd3=diff(qd(:,3));
tn=linspace(0,Tf,(kmax-1));

subplot(2,2,1);
plot(tn,dif_qd1);
ylabel('dqd1/dt (rad)');
xlabel('time t (sec)');

subplot(2,2,2);
plot(tn,dif_qd2);
ylabel('dqd2/dt (rad)');
xlabel('time t (sec)');

subplot(2,2,3);
plot(tn,dif_qd3);
ylabel('dqd3/dt (rad)');
xlabel('time t (sec)'); 
%% robotic manipulator trajectory diagram
fig5 = figure;
axis([-30 30 -30 30 -40 40]) %% limits of 3 axis
grid on
axis on
hold on
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
plot3(xd,yd,zd,'rs'); 
dtk=1000;  %%plotting per dtk samples for the simulation
plot3([0],[0],[0],'o'); 
for tk=1:dtk:kmax,
pause(0.1); %% pausing
 plot3([0,xd00(tk)],[0,yd00(tk)],[0,zd00(tk)]);
 plot3([xd00(tk)],[yd00(tk)],[zd00(tk)],'o');
 plot3([xd00(tk),xd01(tk)],[yd00(tk),yd01(tk)],[zd00(tk),zd01(tk)]);
 plot3([xd01(tk)],[yd01(tk)],[zd01(tk)],'o');
 plot3([xd01(tk),xd02(tk)],[yd01(tk),yd02(tk)],[zd01(tk),zd02(tk)]);
 plot3([xd02(tk)],[yd02(tk)],[zd02(tk)],'o');
 plot3([xd02(tk),xd03(tk)],[yd02(tk),yd03(tk)],[zd02(tk),zd03(tk)]);
 plot3([xd03(tk)],[yd03(tk)],[zd03(tk)],'y*');
 plot3([xd03(tk),xd0E(tk)],[yd03(tk),yd0E(tk)],[zd03(tk),zd0E(tk)]);
 plot3([xd0E(tk)],[yd0E(tk)],[zd0E(tk)],'g+');
end 