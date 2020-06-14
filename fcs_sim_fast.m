%% 'fcs_sim_fast.m' is a script that initiates the efficient simulation of fcs experiment
% In particular, smaller time steps are taken near the confocal volume,
% and larger time step for other region.

% 'fcs_sim_fast.m' writes the intensity trajectories into a file
% 'intensity.bin' 

%% fcs efficient simulation driver

ns = 100;
fprintf("-Number of simulations: %d\n", ns);

dt = 1e-06;
fprintf("Time step dt: %g\n", dt);

T = 1;
fprintf("T: %g\n", T);
taumax = 2e-04;
ntau = floor(taumax/dt);

D = 1.3e-10;
fprintf("-Diffusivity: %g\n", D);

r = 2e-7;
fprintf("Optical system r: %g\n", r);
l = 1e-6;
fprintf("Optical system l: %g\n", l);

br = 8*r;
fprintf("Box length: %g\n", 2*br);
bl = 8*l;
fprintf("Box height: %g\n", 2*bl);
conc = 0.02;
fprintf("Concentration: %g\n", conc);

% derived simulation parameters 
vol_sol = 8*br*br*bl;
avogadro = 6.0221409e+23;
npart_d = conc*avogadro*(1e-06)*1000*vol_sol;
npart_i = round(npart_d);
fprintf("Number of particles %d\n", npart_i);

vol_det = pi * sqrt(pi) * r*r*l;
neff = conc*avogadro*(1e-06)*1000*vol_det;
fprintf("-Effective number of particle in illuminatd region: %g\n", neff);

ntime = round((T+taumax)/dt); 
fprintf("-The number of time points: %d\n", ntime);
nt = ntime - ntau;

intensity_eff = zeros(ntime, ns);

%% efficient simulation
a = 3*r;
c = 3*l;
m = 20;
vol_E = pi*a*a*c*4/3;
ratio = vol_E/vol_sol;

fprintf('--------------------efficient simulation---------------\n')
tic
for j = 1:ns
    %intensity_eff(:,j) =  get_intensity_ellip(dt, ntime, D, br, bl, npart_i, r, l, a, c, m);
    %intensity_eff(:,j) =  get_intensity_fast(dt, ntime, D, br, bl, npart_i, r, l);
    intensity_eff(:,j) = get_intensity_eff(dt, ntime, D, br, bl, ... 
		npart_i, r, l, m);
end
toc

t_dat = (1:ntime)*dt;

fileID = fopen('intensity.bin','w');
fwrite(fileID,intensity_eff,'double');
fclose(fileID);

clear intensity_eff
%%
%{
fileID = fopen('intensity.bin');
inten_vec = fread(fileID, ns*ntime, 'double');
fclose(fileID);
%%
intensity_recover = reshape(inten_vec, [ntime, ns]);

%% Correlation
Corr_dat = zeros(ns, ntau);
Intensity_1 = intensity_eff(1:ns, 1:nt); % ns by nt matrix
for k = 1:ntau
    Intensity_k = intensity_eff(1:ns, k:(k+nt-1));
    % dot(A, B, 2) computes the inner product of the rows of matrix A, B
    % where A, B are m by n matrices, dot(A, B, 2) is a m by 1 vector
    Corr_dat(:,k) = dot(Intensity_1, Intensity_k, 2)/ntime;
end
%}