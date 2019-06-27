l = 5.6983;

openfemm
opendocument('simulacija.FEM')

%mi_selectgroup(2)
%mi_probdef(2.5, 'centimeters', 'planar', 1E-8, l, 30, 0);
%mi_analyze(1);
mi_loadsolution;
mo_seteditmode('area');
mo_selectblock(0, 3.5); 

for k = 1:1:48
    mo_selectblock(4.3*cos(angle(k)), 4.3*sin(angle(k)));    
end

M = mo_blockintegral(22)