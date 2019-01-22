%% Basic symmetry analysis
% We will demonstrate the symmetry analysis capabilities of SpinW in this
% script. We determine the allowed elements of the Dzyaloshinskii-Moriya
% vector on a bond that is perpendicular to a mirror plane. The result is
% well known and can be compared to the SpinW output. We use here
% explicitly the mirror plane as the only symmetry operator. Also we will
% draw the plane to demonstrate the new plotting functions of SpinW.

model = spinw;
model.genlattice('lat_const',[6 4 4],'sym','-x,y,z','label','m_x')
model.addatom('r',[3/8 1/2 1/2],'S',1)
plot(model)

% draw mirror plane
Rlu = permute([1/2 0 0;1/2 1 0;1/2 1 1;1/2 0 1],[2 3 1]);
swplot.plot('type','polyhedron','position',Rlu,'color','pink','alpha',0.5)

model.gencoupling
model.addmatrix('label','DM1','value',1);
model.addcoupling('mat','DM1','bond',1)
plot(model)
model.getmatrix('mat','DM1');

% Look at the figure from a special angle using the swplot.transform
% function and recenter the figure using swplot.translate.
R = sw_rotmatd([1 1 0],45);
swplot.transform(R)
swplot.translate
    
%% Draw the 2 allowed DM components
% We add two DM vector that span the allowed plane.

model.addmatrix('label','DM1','value',[0 1 0],'color','blue')
model.addmatrix('label','DM2','value',[0 0 1],'color','purple')
model.addcoupling('mat','DM1','bond',1)
model.addcoupling('mat','DM2','bond',1)
plot(model,'bondScale',1)

%% Mirror plane including the bond
% The second example show the case when the mirror plane includes the bond.

model = spinw;
model.genlattice('lat_const',[4 6 4],'sym','-x,y,z','label','m_x')
model.addatom('r',[1/2 3/8 1/2],'S',1)
model.addatom('r',[1/2 5/8 1/2],'S',1)
plot(model)

% Lets draw mirror plane.
Rlu = permute([1/2 0 0;1/2 1 0;1/2 1 1;1/2 0 1],[2 3 1]);
swplot.plot('type','polyhedron','position',cat(3,Rlu),'color','pink',...
    'alpha',0.5,'replace',true)

model.quickham(1)
plot(model)
model.getmatrix('mat','J1');

R = sw_rotmatd([0.5 1 0],45);
swplot.transform(R)
swplot.translate

%% DM interaction using P4 space group
% SpinW not only determines the symmetry allowed matrix elements of the
% Hamiltonian, but it also transforms the exchange (anisotropy and g)
% matrices using the symmetry operators. This example shows this.

cryst=spinw;
cryst.genlattice('sym','P 4','lat_const',[8 8 6])
cryst.addatom('r',[1/4 1/4 0],'S',1);
cryst.gencoupling;
cryst.addmatrix('label','D','value',[1 -1 0])
cryst.addcoupling('mat','D','bond',2)

plot(cryst,'range',[1 1 1/2])

%% DM interaction using no space group (P0)

cryst=spinw;
cryst.genlattice('lat_const',[8 8 6])
cryst.addatom('r',[1/4 1/4 0],'S',1);
cryst.addatom('r',[3/4 1/4 0],'S',1);
cryst.addatom('r',[1/4 3/4 0],'S',1);
cryst.addatom('r',[3/4 3/4 0],'S',1);
cryst.gencoupling
cryst.addmatrix('label','D','value',[1 -1 0])
cryst.addcoupling('mat','D','bond',1,'subIdx',[3 5 7:8])

plot(cryst,'range',[1 1 1/2],'atomLegend',false)

%% single ion anisotropy using P4 space group

cryst=spinw;
cryst.genlattice('sym','P 4','lat_const',[8 8 6])
cryst.addatom('r',[1/4 1/4 0],'S',1)
cryst.gencoupling
cryst.addmatrix('label','A','value',1-eye(3))
cryst.addaniso('A')

plot(cryst)

%% single ion anisotropy using no space group

cryst=spinw;
cryst.genlattice('lat_const',[8 8 6])
cryst.addatom('r',[1/4 1/4 0],'S',1)
cryst.addatom('r',[3/4 1/4 0],'S',1)
cryst.addatom('r',[1/4 3/4 0],'S',1)
cryst.addatom('r',[3/4 3/4 0],'S',1)
cryst.gencoupling
cryst.addmatrix('label','A','value',1-eye(3))
cryst.addaniso('A')

plot(cryst,'atomLegend',false)


%% Spin wave spectrum of the honeycomb lattice Na2IrO3
% We reproduce structure presented in the paper
% S. K. Choi, et al. PRL, 108(12), 127204 (2012),
% [[http://link.aps.org/doi/10.1103/PhysRevLett.108.127204]].


%% Na2IrO3 crystal structure
% We define the crystal structure using the space group C2/m, and taking
% the crystallographic parameters at 300 K (parameters are only slightly
% different at 5 K) and we add not only the magnetic Ir4+ ions with effective
% spin quantum number of 1/2 but also the non-magnetic atoms for plotting
% the structure.

nairo = spinw;
nairo.genlattice('lat_const',[5.427 9.395 5.614],'angled',[90 109.037 90],'spgr','C 2/m')
nairo.addatom('label','MIr4','r',[1/2; 0.167; 0],'S',1/2,'color','DarkCyan');
nairo.addatom('r',[0 1/2 1/2;0 0 0.340; 0 1/2 1/2],'S',[0 0 0],'label',{'Na1 Na1+' 'Na2 Na1+' 'Na3 Na1+'},'color',{'lightGray' 'lightGray' 'lightGray'});
nairo.addatom('r',[0.748 0.711; 0.178 0; 0.789 0.204],'S',[0 0],'label',{'O1 O2-', 'O2 O2-'},'color',{'r' 'r'});
plot(nairo,'baseShift',[-2;0;0])
swplot.zoom(1.3)

% We generate all bonds up to 8 Angstrom length.
nairo.gencoupling('maxDistance',8);

%% Magnetic Hamiltonian
% The Kitaev model is a spin model on the honeycomb lattice. It assumes
% that there are 3 type of exchange interaction on the three different
% types of first neighbor bonds running along the three different
% directions within the lattice plane. Exchange couplings
K   = 1;
Jxx = diag([K 0 0]);
Jyy = diag([0 K 0]);
% and
Jzz = diag([0 0 K]);
% are assigned to the three types of bonds.
%
% Let's implement this model in SpinW for Na2IrO3!
%
% We define the three anisotropic exchange interactions (Jxx, Jyy and Jzz)
% and assign them according to the paper.

nairo.addmatrix('label','Jxx','value',Jxx,'color','r');
nairo.addmatrix('label','Jyy','value',Jyy,'color','g');
nairo.addmatrix('label','Jzz','value',Jzz,'color','b');

% We start assigning the Kitaev term Jxx to the first neighbor bonds.
nairo.addcoupling('mat','Jxx','bond',1);

plot(nairo,'range',[2 2 0.5],'atomMode','mag','cellMode','inside',...
    'atomLegend',false,'cellcolor','gray','bondMode','line','bondLinewidth0',2)

%% Do you see any problem here?
%
% Let's do symmetry analysis
nairo.getmatrix('mat','Jxx')
% Looks fine, but what to do???
%
% Remove previous bond-matrix assignments:
nairo.gencoupling
% Let's break the lattice symmetry:
nairo.addcoupling('mat','Jxx','bond',1,'subidx',[3 4]);
nairo.addcoupling('mat','Jyy','bond',1,'subidx',[1 2]);

plot(nairo,'range',[2 2 0.5],'atomMode','mag','cellMode','inside',...
    'atomLegend',false,'cellcolor','gray','bondMode','line','bondLinewidth0',2)
% Also add Jzz
nairo.addcoupling('mat','Jzz','bond',2);
% Plot Kitaev couplings only.
plot(nairo,'range',[2 2 0.5],'atomMode','mag','cellMode','inside',...
    'atomLegend',false,'cellcolor','gray','bondMode','line','bondLinewidth0',2)
swplot.zoom(1.4)

%% Spin wave dispersion
% This will be possible after the magnetic structure presentation. If you
% feel adventurous try to generate the magnetic structure and plot the
% dispersion over the Q ranges in the paper.





