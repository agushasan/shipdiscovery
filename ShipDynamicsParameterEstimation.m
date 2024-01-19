%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 10;
dt  = 0.001;
t   = dt:dt:tf;

nn = 3;
mm = 27;

%% system description
A = eye(6);
B_t = [1 0;0 0; 0 1];
M = [55 0 0;0 66 12;0 12 5];
B = [0 0;0 0;0 0;inv(M)*B_t];
C = [0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];
Cb = eye(3);

%% noise
QF = 0.001*eye(rank(A));
RF = 0.001*eye(rank(C));

%% state initialization
x        = [0;0;0;0;0;0];
y        = [0;0;0];
vbar     = [0;0;0];
thetabar = zeros(mm,1);
 
%% known paramaters
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);
m   = m22*m33-m23*m32;

%% True parameters
Xuu = -4.2;
Yvv = -9.3;
Nrr = -0.5;
Xu  = -15.6;
Yv  = -25.1;
Yr  = -0.6;
Nv  = -1.2;
Nr  = -2.5;

%% initial control inputs
u     = [100 1]';

%% for plotting
uArray          = [];
xArray          = [];
yArray          = [];
vbarArray       = [];
thetabarArray   = [];

%% Initialization for estimator

lambdav = 0.995;
lambdat = 0.999;
Rv      = 0.001*eye(nn);
Rt      = 0.001*eye(nn);
Pv      = 100000*eye(nn);
Pt      = 100000*eye(mm);
Gamma   = zeros(nn,mm);

%% simulation
for i=1:(tf/dt)
    
    u     = [100*sin(i*dt) 1*cos(i*dt)]';

    uArray         = [uArray u];
    xArray         = [xArray x];
    yArray         = [yArray y];
    vbarArray      = [vbarArray vbar];
    thetabarArray  = [thetabarArray thetabar]; 

    c13 = -m22*x(5)-((m23+m32)/2)*x(6);
    c23 = m11*x(4);
    Cv  = [0 0 c13; 0 0 c23;-c13 -c23 0];
    Dv  = -[Xu+Xuu*abs(x(4)) 0 0;0 Yv+Yvv*abs(x(5)) Yr;0 Nv Nr+Nrr*abs(x(6))];
    
    x = A*x+dt*[cos(x(3))*x(4)-sin(x(3))*x(5);sin(x(3))*x(4)+cos(x(3))*x(5);x(6);-inv(M)*(Cv+Dv)*[x(4);x(5);x(6)]]+B*dt*u;%+QF*dt*randn(6,1);
    y = C*x;%+RF*dt*randn(rank(C),1);

%     Phi = [y(1) abs(y(1))*y(1) 0 0 0 0 0 0 0 0;
%            0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3) 0 0 0 0;
%            0 0 0 0 0 0 y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3)];
%     Phi = [y(1) y(2) y(3) abs(y(1))*y(1) abs(y(2))*y(2) abs(y(3))*y(3) zeros(1,12);
%            zeros(1,6) y(1) y(2) y(3) abs(y(1))*y(1) abs(y(2))*y(2) abs(y(3))*y(3) zeros(1,6);
%            zeros(1,12) y(1) y(2) y(3) abs(y(1))*y(1) abs(y(2))*y(2) abs(y(3))*y(3)];

    Phi = [y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3 zeros(1,18);
           zeros(1,9) y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3 zeros(1,9);
           zeros(1,18) y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)^3 y(2)^3 y(3)^3];

    bk  = dt*inv(M)*(B_t*u-Cv*y);
    
    % Estimation using adaptive observer
    Kv = Pv*Cb'*inv(Cb*Pv*Cb'+Rv);
    Kt = Pt*Gamma'*Cb'*inv(Cb*Gamma*Pt*Gamma'*Cb'+Rt);
    Gamma = (eye(nn)-Kv*Cb)*Gamma;

    vbar = vbar+(Kv+Gamma*Kt)*(y-Cb*vbar);
    thetabar = thetabar-Kt*(y-Cb*vbar);

    vbar = eye(nn)*vbar+bk+Phi*thetabar;

    thetabar = thetabar;
    Pv = (1/lambdav)*eye(nn)*(eye(nn)-Kv*Cb)*Pv*eye(nn);
    Pt = (1/lambdat)*(eye(mm)-Kt*Cb*Gamma)*Pt;
    Gamma = eye(rank(Cb))*Gamma-Phi;

end

Temp1 = inv([-dt*m23/m dt*m33/m;dt*m22/m -dt*m32/m])*[thetabarArray(3,:);thetabarArray(7,:)];
Temp2 = inv([-dt*m23/m dt*m33/m;dt*m22/m -dt*m32/m])*[thetabarArray(4,:);thetabarArray(8,:)];

figure(1)
yyaxis left
plot(t,uArray(1,:), ':b', 'LineWidth', 16)
yyaxis right
plot(t,uArray(2,:), ':r', 'LineWidth', 16)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
legend('\tau_u','\tau_r','FontSize',72)
xlabel('time (s)','FontSize',72)

figure(2)
subplot(3,1,1)
plot(t,yArray(1,:), '-k', 'LineWidth', 16)
hold on;
plot(t,vbarArray(1,:), ':g', 'LineWidth', 16)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('u [m/s]','FontSize',72)
subplot(3,1,2)
plot(t,yArray(2,:), '-k', 'LineWidth', 16)
hold on;
plot(t,vbarArray(2,:), ':g', 'LineWidth', 16)
grid on;
grid minor;
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylabel('v [m/s]','FontSize',72)
legend('measurement','estimated','FontSize',72)
subplot(3,1,3)
plot(t,yArray(3,:), '-k', 'LineWidth', 16)
hold on;
plot(t,vbarArray(3,:), ':g', 'LineWidth', 16)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('r [rad/s]','FontSize',72)
xlabel('time (s)','FontSize',72)

figure(3)
subplot(4,2,1)
plot(t,Xu*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,m11*thetabarArray(1,:)/dt, ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Xu-2 Xu+2]);
ylabel('X_u','FontSize',72)
subplot(4,2,2)
plot(t,Xuu*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,m11*thetabarArray(2,:)/dt, ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('X_{uu}','FontSize',72)
legend('true parameter','estimated parameter','FontSize',48)
ylim([Xuu-2 Xuu+2]);
subplot(4,2,3)
plot(t,Nv*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,Temp1(1,:), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nv-2 Nv+2]);
ylabel('N_v','FontSize',72)
subplot(4,2,4)
plot(t,Yv*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,Temp1(2,:), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_v','FontSize',72)
ylim([Yv-2 Yv+2]);
subplot(4,2,5)
plot(t,Nr*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,Temp2(1,:), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([Nr-2 Nr+2]);
ylabel('N_r','FontSize',72)
subplot(4,2,6)
plot(t,Yr*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,Temp2(2,:), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('Y_r','FontSize',72)
ylim([Yr-2 Yr+2]);
subplot(4,2,7)
plot(t,Yvv*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,m*thetabarArray(5,:)/(dt*m33), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('time (s)','FontSize',72)
ylim([Yvv-2 Yvv+2]);
ylabel('Y_{vv}','FontSize',72)
subplot(4,2,8)
plot(t,Nrr*ones(length(t),1), '-k', 'LineWidth', 14)
hold on;
plot(t,-m*thetabarArray(6,:)/(dt*m23), ':g', 'LineWidth', 14)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('N_{rr}','FontSize',72)
xlabel('time (s)','FontSize',72)
ylim([Nrr-2 Nrr+2]);