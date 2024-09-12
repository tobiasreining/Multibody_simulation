    clear all
close all
set(0,'defaultFigureWindowStyle','default')

l=10; 
m=15;
mu_min=(-3^0.5+2)/2;
mu=0.95*mu_min;
g=9.81;

%Definition der Variablen und Ableitungen
syms phi(t) t 
Dphi=diff(phi,t,1);
Dphi2=diff(phi,t,2);

%Bewegungsgleichung
eq_kontakt=Dphi2*(mu*l/4*sin(2*phi)-l/3)-mu*l/2*sin(phi)^2*Dphi^2+g*(mu*sin(phi)-cos(phi)/2)==0;

%Bewegungsgleichung lösen
[eqs,vars]=reduceDifferentialOrder(eq_kontakt,phi(t));
[M,F]=massMatrixForm(eqs,vars);
f=M\F;
fh=odeFunction(f,vars);
phi0=[5*pi/12,0];  %Anfangsbedingung
opt=odeset('maxstep',0.01);
tspan=[0 5]; %Zeitbereich für Integration
solver=@ode45;
sol=solver(@(t,phi) fh(t,phi),tspan,phi0,opt);

%Daten lesen
time=sol.x;
angle=sol.y;
p=angle(1,:);
dp=angle(2,:);

%Gleichungen für Größen aus DS und SPS, bei Kontakt
ddp=(mu*sin(p).^2.*dp.^2-g*(mu*sin(p)-cos(p)/2))./(mu*l/4*sin(2*p)-1/3*l);
FSide=m*mu*(g+l/2*(cos(p).*ddp-sin(p).*dp.^2))-1/2*m*l*(sin(p).*ddp+cos(p).*dp.^2);
R=m*mu*(g+l/2*(cos(p).*ddp-sin(p).*dp.^2));
N=m*(g+l/2*(cos(p).*ddp-sin(p).*dp.^2));

%Werte Bei Kontaktverlust
loc_lose=find(abs(FSide)==min(abs(FSide))); %Ort in Array bei Kontaktverlust
t_lose=time(loc_lose);  %findet Zeitpunkt des Kontaktverlust
p_lose=p(loc_lose);
dp_lose=dp(loc_lose);

%Begrenzen der Werte bis Kontaktverlust
time_kontakt=time(1:loc_lose);
p_kontakt=p(1:loc_lose);
dp_kontakt=dp(1:loc_lose);
ddp_kontakt=ddp(1:loc_lose);
FSide_kontakt=FSide(1:loc_lose);
R_kontakt=R(1:loc_lose);
N_kontakt=N(1:loc_lose);

%Neue Bewegungsgleichung nach Kontaktverlust
eq_lose=Dphi2*(mu*l/4*sin(2*phi)-cos(phi)^2*l/2-l/6)+Dphi^2*(sin(2*phi)*l/4-sin(phi)^2*mu*l/2)+g*(mu*sin(phi)-cos(phi))==0;

%Lösen der neuen Bewegungsgleichung
[eqs,vars]=reduceDifferentialOrder(eq_lose,phi(t));
[M,F]=massMatrixForm(eqs,vars);
f=M\F;
fh=odeFunction(f,vars);
phi0=[p_lose,dp_lose]; %AB für DGL2 sind Werte von DGL1 bei Kontaktverlust
tspan=[0 1]; 
sol=solver(@(t,phi) fh(t,phi),tspan,phi0,opt);

%Neue Daten lesen
time_lose=sol.x(2:end);
time_lose=time_lose + t_lose; %Damit das neue 'time' nicht bei 0 beginnt
angle_lose=sol.y;
p_lose=angle_lose(1,2:end);
dp_lose=angle_lose(2,2:end);

%Gleichungen für Kräfte, aus DS und SPS, nach Kontaktverlust
ddp_lose=(dp_lose.^2.*(sin(2*p_lose)*l/4-sin(p_lose).^2*mu*l/2)+g*(mu*sin(p_lose)-cos(p_lose)))./(-mu*l/4*sin(2*p_lose)+cos(p_lose).^2*l/2+l/6);
R_lose=m*mu*(g+l/2*(cos(p_lose).*ddp_lose-sin(p_lose).*dp_lose.^2));
N_lose=m*(g+l/2*(cos(p_lose).*ddp_lose-sin(p_lose).*dp_lose.^2));

%Aneinanderhängen der Daten
time_gesamt=[time_kontakt time_lose];
p_gesamt=[p_kontakt p_lose];
dp_gesamt=[dp_kontakt dp_lose];
FSide_gesamt=[FSide_kontakt zeros(1, length(time_lose))]; %F ist 0 nach Kontaktverlust
loc_end=find(abs(p_gesamt)==min(abs(p_gesamt))); %array location am Ende
t_gekuerzt=time_gesamt(1:loc_end); %gekürzter Zeitintervall für R und N aus Darstellungsgründen
p_gekuerzt=p_gesamt(1:loc_end);
R_gesamt=[R_kontakt R_lose];
R_gesamt=R_gesamt(1:loc_end);
N_gesamt=[N_kontakt N_lose];
N_gesamt=N_gesamt(1:loc_end);

%SIMPACK Vergleichsdaten
Simpack=load('Leiter095.mat');
F_Simpack=Simpack.timeInt.forceForce.SF_side.abs.values;
R_Simpack=Simpack.timeInt.forceForce.SF_base_friction.abs.values;
N_Simpack=Simpack.timeInt.constrForce.SL_base.la_001.values;
alpha_Simpack=Simpack.timeInt.bodyPosRot.SB_Body1.alpha.values;
t_Simpack=Simpack.timeInt.time.values;

%Plotten
figure(1)
subplot(2,2,1)
plot(t_Simpack,alpha_Simpack)
title('\alpha(t) für \mu = 0.95* \mu_hmin')
hold on
plot(time_gesamt,p_gesamt)
ylim([0,1.5])
legend('SIMPACK','MATLAB')

subplot(2,2,2)
plot(t_Simpack,F_Simpack)
title('Kraft an der Wand über t für \mu = 0.95* \mu_hmin')
hold on
plot(time_gesamt,FSide_gesamt)
ylim([0,50])
legend('SIMPACK','MATLAB')

subplot(2,2,3)
plot(t_Simpack,R_Simpack)
title('Reibkraft über t für \mu = 0.95* \mu_hmin')
hold on
plot(t_gekuerzt,R_gesamt)
ylim([0,50])
legend('SIMPACK','MATLAB')

subplot(2,2,4)
plot(t_Simpack,N_Simpack)
title('Normalkraft über t für \mu = 0.95* \mu_hmin')
hold on
plot(t_gekuerzt,N_gesamt)
ylim([0,160])
legend('SIMPACK','MATLAB')

figure(2)
subplot(2,2,1)
plot(alpha_Simpack,F_Simpack)
title('Kraft an der Wand über \alpha für \mu = 0.95* \mu_hmin')
hold on
plot(p_gesamt,FSide_gesamt)
xlim([0,1.4])
ylim([0,50])
legend('SIMPACK','MATLAB')

subplot(2,2,2)
plot(alpha_Simpack,R_Simpack)
title('Reibkraft über \alpha für \mu = 0.95* \mu_hmin')
hold on
plot(p_gekuerzt,R_gesamt)
xlim([0,1.4])
legend('SIMPACK','MATLAB')

subplot(2,2,3)
plot(alpha_Simpack,N_Simpack)
title('Normalkraft über \alpha für \mu = 0.95* \mu_hmin')
hold on
plot(p_gekuerzt,N_gesamt)
xlim([0,1.4])
legend('SIMPACK','MATLAB')

close all %Damit im PDF keine Bilder sind