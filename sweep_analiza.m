%% Podaci iz proracuna za projekat
Dse0 = 16.3877; 
Dsyi = 14.9103;
Dryi= 0.0815;
Dri= 0.0657;
hys0 = 0.7387;
l = 5.6983;
Roal = 2720;
Rocu = 8933;
Rofe = 7250;
S1prs = 0.65e-6;
S1prr = 2.4387e-5;
Qr = 48;
Qs = 36;
wts = 0.0034;
wtr = 0.0024;
hts = 0.0271;
htr = 0.0067;
kispune = 0.55;
%%
hys = (0.8 : 0.02 : 1.2) * hys0;
Dse = Dsyi + 2*hys;
f = (0.01 : 0.01 : 0.1)  * 50;
angle = (7.5*(1:1:Qr) - 3.6)*pi/180;
%M = zeros(length(Dse), length(f));
M = zeros(1, length(f));
%Masau=zeros(length(Dse),length(f));
Masau = zeros(1, length(f));
Masar = ((Dryi^2-Dri^2)*pi/4+wtr*htr*Qr)*l*Rofe+S1prr*l*Roal*Qr;

%%
openfemm
opendocument('simulacija.FEM')

for i = 1 : 1 : length(Dse)
    mi_selectgroup(2)
    mi_scale2(0, 0, Dse(i)/Dse0, 4);
    
    for j = 1 : 1 : length(f)
        mi_probdef(f(j), 'centimeters', 'planar', 1E-8, l, 30, 0);
        mi_analyze(0);
        
        mi_loadsolution;
        mo_seteditmode('area');
        mo_selectblock(0, 3.5);  
        for k = 1:1:Qr
             mo_selectblock(4.3*cos(angle(k)), 4.3*sin(angle(k)));    
        end
        
        Masas=((Dse(i)^2-Dsyi^2)*pi/4+wts*hts*Qs)*l*Rofe+kispune*S1prs*l*Rocu*Qs;
        Masau(i,j)=Masar+Masas;
        M(1, j) = mo_blockintegral(22)
    end
    
    mi_selectgroup(2);
    mi_scale2(0, 0, Dse0/Dse(i), 4);
end