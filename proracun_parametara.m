close all, clear all, clc

%proracun osnovnih dimenzija rotora
SIGMAftan=13990;
Pn=1100;
In=3.25;
n=930;
p=3;
q=3;
Bm=0.7;
Un=400;
omegan=2*pi*n/60;
Vr=Pn/(2*SIGMAftan*omegan);
X=0.6; %pi*p^(1/3)/(2*p);
D=(4*Vr/(X*pi))^(1/3);
l=X*D;
%izbor broja zlebova i vrste namotaja 
Qs=36;
taus=pi*D/Qs;
taup=pi*D/(2*p);
%proracun sirine medjugvozdja
delta=0.8e-3; %(0.18+0.006*Pn^0.4)/1000;
%proracun potrebnog broja navojaka faznog namotaja 
Ufn=Un/sqrt(3);
E=0.75*Ufn;
alfai=2/pi;
f=50;
w=2*pi*f;
z=Qs/(2*p);
m=Qs/(2*p*q);
y=5;
kw1=sin(m*pi/(2*z))/(m*sin(pi/(2*z)))*sin(y*taus*pi/(taup*2));
Ns=sqrt(2)*E/(w*kw1*alfai*taup*l*Bm);
as=1;
Ns=696; %ceil(Ns);
Zqs=2*as*q*Ns/Qs;
%provera poduzne gustine struje
t=0.69;
A=sqrt(2)*SIGMAftan/(Bm*t);
A2=2*In*Ns*q/(pi*D);
j=(A-A2)/A2*100
%proracun dimenzija zubaca i zlebova 
Qr=48;
taur=pi*D/Qr;
Bts=1.7;
Btr=1.8;
wts= Bm*taus/Bts; % Bm*taus/Bts+1/2*Bm*taus/Bts;
wtr= Bm*taur/Btr; % Bm*taur/Btr+1/2*Bm*taur/Btr;
wss=taus-wts;
wsr=taur-wtr;
Kss=0.6;
Ksr=1;
Jmax=8000000;
S1prr=Zqs*Qs*In*t/(as*Qr*Jmax);
hsr=Zqs*Qs*In*t/(as*Qr*Jmax*wsr);
S1prs=0.65e-6;
hss=Zqs*S1prs/(wss*Kss);
Sqs=Kss*wss*hss;
Sqr=Ksr*wsr*hsr;
hs=0.05*hss;
hr=0.05*hsr;
hts=hss+hs;
htr=hsr+hr;
%provera opravdanosti zanemarenja zasicenja u magnetskom kolu 
mi0=4*pi*10^(-7);
Hts=278;
Htr=278;
Utr=Htr*htr;
Uts=Hts*hts;
Udelta=Bm*delta/mi0;
Ksat=(Uts+Utr)/Udelta;
%proracun sirine jarma statora i rotora 
F=alfai*Bm*taup*l;
Bys = 1.5;
Byr = 1.4;
hys=F/(2*Bys*l);
hyr=F/(2*Byr*l);
Dryi=D-2*htr;
Dri=Dryi-2*hyr;
Dsyi=D+2*hts;
Dse=Dsyi+2*hys;
%proracun gubitaka u magnetskom kolu masine 
pfeys=2.42;
pfets=2.2;
gamafe=7650;
Pfes=(pfeys*pi*(Dse^2-Dsyi^2)/4+pfets*Qs*hts*wts)*l*gamafe;
%proracun gubitaka u namotajima 
ktetas=1.4;
sigmacu=57000000;
Pcus=3*ktetas*Ns*(2*l+3.8*y*taus)*In^2/(S1prs*sigmacu)
ktetar=1.5;
sigmaal=35000000;
alfar=2*pi*p/Qr;
Sring=Sqr/(2*sin(alfar/2));
Rr=Qr/3*(ktetar*l/(Sqr*sigmaal)+ktetar*taur/(2*Sring*sigmaal*(sin(alfar/2))^2));
Pcur=3*Rr*(Zqs*Qs*In*t/(as*Qr))^2
%proracun struje magnecenja
bs=0.4*wss;
br=0.4*wsr;
kc=taus/(taus-(bs/delta)/(5+bs/delta)*bs);
deltaef=kc*delta;
Im=pi*p*Bm*deltaef/(mi0*q*kw1*Ns*sqrt(2));
%proracun parametara masine
krs=4*q*(kw1*Ns)^2/Qr;
Rr1=krs*Rr;
Lm=6*mi0*taup*l*(kw1*Ns)^2/(pi^2*p*deltaef);
Rfe=3*Ufn^2/Pfes;
k1=1-9/16*(1-y*taus/taup);
k2=1-3/4*(1-y*taus/taup);
Lqs=4*q*mi0*l*Ns^2*(k1*hss/(3*wss)+k2*hs/wss)/Qs;
Lzs=4*q*mi0*l*Ns^2*k2*(5*delta/bs)/(5+4*delta/bs)/Qs;
Lws=4*q*m*Ns^2*mi0*2*y*taus*0.3/Qs;
koef_sigma=0;
for i=3:2:53,
    kwi = sin(i*m*pi/(2*z))/(m*sin(i*pi/(2*z)))*sin(i*y*taus*pi/(taup*2)); 
    koef_sigma=koef_sigma+kwi^2/i^2/kw1^2;
end
sigmas=koef_sigma;
Lsigmas=sigmas*Lm+Lqs+Lzs+Lws;
Ldeltar1=pi^2/3*(p/Qr)^2*Lm;
Lqr1=krs*mi0*l*(hsr/(3*wsr)+hr/wsr);
Lzr1=krs*mi0*l*(5*delta/br)/(5+4*delta/br);
Lwr1=krs*mi0*Qr*(0.18*pi*D/(2*p))/(q*p^2*3);
Lsigmar1=Ldeltar1+Lqr1+Lzr1+Lwr1;
Xm=w*Lm;
Xsigmas=w*Lsigmas;
Xsigmar = w*Lsigmar1;
Xdeltar1=w*Ldeltar1;
Rs = Pcus/3/In^2;
Rr = Pcur/3/In^2;

% Parametri u relativnim jedinicama:
Zb = Ufn/In;
xm_pu = Xm/Zb
xsigmas_pu = Xsigmas/Zb
xsigmar_pu = Xsigmar/Zb
rs_pu = Rs/Zb
rr_pu = Rr/Zb

%% Slicica:
angle = 0:0.001:2*pi;
figure, plot(Dri/2*cos(angle),Dri/2*sin(angle)), hold on
plot(Dryi/2*cos(angle),Dryi/2*sin(angle),':')
plot(Dsyi/2*cos(angle),Dsyi/2*sin(angle),':')
plot(Dse/2*cos(angle),Dse/2*sin(angle))
plot(D/2*cos(angle),D/2*sin(angle))
plot((D/2+delta)*cos(angle),(D/2+delta)*sin(angle))
ss_angle = wss/taus*2*pi/Qs;
% statorski zlebovi:
for i = 1:Qs
    plot([(D/2+delta)*cos((i-1)/Qs*2*pi+ss_angle/2),(D/2+delta+hss)*cos((i-1)/Qs*2*pi+ss_angle/2)],...
        [(D/2+delta)*sin((i-1)/Qs*2*pi+ss_angle/2),(D/2+delta+hss)*sin((i-1)/Qs*2*pi+ss_angle/2)])
    plot([(D/2+delta)*cos((i-1)/Qs*2*pi-ss_angle/2),(D/2+delta+hss)*cos((i-1)/Qs*2*pi-ss_angle/2)],...
        [(D/2+delta)*sin((i-1)/Qs*2*pi-ss_angle/2),(D/2+delta+hss)*sin((i-1)/Qs*2*pi-ss_angle/2)])
    plot([(D/2+delta+hss)*cos((i-1)/Qs*2*pi+ss_angle/2),(D/2+delta+hss)*cos((i-1)/Qs*2*pi-ss_angle/2)],...
        [(D/2+delta+hss)*sin((i-1)/Qs*2*pi+ss_angle/2),(D/2+delta+hss)*sin((i-1)/Qs*2*pi-ss_angle/2)])
end
% rotorski zlebovi:
sr_angle = wsr/taur*2*pi/Qr;
for i = 1:Qr
    plot([D/2*cos((i-1)/Qr*2*pi+sr_angle/2),(D/2-hsr)*cos((i-1)/Qr*2*pi+sr_angle/2)],...
        [D/2*sin((i-1)/Qr*2*pi+sr_angle/2),(D/2-hsr)*sin((i-1)/Qr*2*pi+sr_angle/2)])
    plot([D/2*cos((i-1)/Qr*2*pi-sr_angle/2),(D/2-hsr)*cos((i-1)/Qr*2*pi-sr_angle/2)],...
        [D/2*sin((i-1)/Qr*2*pi-sr_angle/2),(D/2-hsr)*sin((i-1)/Qr*2*pi-sr_angle/2)])
    plot([(D/2-hsr)*cos((i-1)/Qr*2*pi+sr_angle/2),(D/2-hsr)*cos((i-1)/Qr*2*pi-sr_angle/2)],...
        [(D/2-hsr)*sin((i-1)/Qr*2*pi+sr_angle/2),(D/2-hsr)*sin((i-1)/Qr*2*pi-sr_angle/2)])
end

xlim([-1.1*Dse/2 1.1*Dse/2])
ylim([-1.1*Dse/2 1.1*Dse/2])
axis equal
