function [intensity_accurate, intensity_efficient] =  get_intensity_coupled_new(dt, nt, D, br, bl, ... 
		npart, r, l, m)
%% efficient simulation of FCS intensity trajectory

% Update the particle position every dt in an important region,
% and the rest of the particles are evolved every (Delta t)
% Delta t = m * dt, m is an integer

%---------------important region-----------------------
% important region: E = {(x1, x2, x3):(x1/a)^2 + (x2/a)^2 + (x3/c)^2 <= 1 }
%  {(x1, x2, x3) : x1^2/a^2 + x2^2/a^2 + x3^2/c^2 <= 1 }

intensity_accurate = zeros(nt, 1);
intensity_efficient = zeros(nt, 1);

sigma_in = sqrt(2*D*dt);
sigma_out = sqrt(m)*sigma_in;

% inside_status: logical array to record if a particle is in the fast region.
%inside_status = zeros(nt, 1);

% wall_status: logical array to record if a particle is near the wall of
% the solution box, and hence need to check reflection.
%wall_status = zeros(nt, 1);

% x_dat(i,:) = (x1, x2, x3) for the ith particle
x_dat = zeros(npart, 3);

% intensity(i) is the intensity from ith particle
%intensity = zeros(nt, 1);

% rho = x1^2/r^2 + x2^2/r^2 + x3^2/l^2
rho_dat = zeros(npart, 1);

a = 3*r;
c = 3*l;

a2 = a*a;
c2 = c*c;
r2 = r*r;
l2 = l*l;

%lmax = floor(nt/m);

% initial position, uniformly distributed
x1_init = 2*br*rand(npart, 1) - br;
x2_init = 2*br*rand(npart, 1) - br;
x3_init = 2*bl*rand(npart, 1) - bl;
x_dat = [x1_init, x2_init, x3_init];

rho_dat = x_dat(:,1).^2/r2 + x_dat(:,2).^2/r2 + x_dat(:,3).^2/l2;
intensity_accurate(1) = sum(exp(-2*rho_dat));


% rho_dat > rho0 iff the particle inside the important region
rho0 = a2/r2;
% rho_dat > rho1 iff we would like to check reflection for that particle
rho1 = (br-r)^2/r2;

% rho_dat > rho0 iff x1^2/a^2 + x2^2 + x3^2 < 1. 
inside_status = rho_dat < rho0;
intensity_efficient(1) = sum(exp(-2*rho_dat(inside_status)));
npart_in = sum(inside_status);

for j = 2:nt
    % update the locations of the particles 
    x_dat = x_dat + sigma_in * randn(npart, 3);
    % check reflection at the wall of the solution box
    ref_ind = x_dat(:, 1) > br;
    x_dat(ref_ind, 1) = 2*br - x_dat(ref_ind,1);
    ref_ind = x_dat(:, 1) < -br;
    x_dat(ref_ind, 1) = -x_dat(ref_ind, 1)-2*br;
    
    ref_ind = x_dat(:, 2) > br;
    x_dat(ref_ind, 2) = 2*br - x_dat(ref_ind, 2);
    ref_ind = x_dat(:, 2) < -br;
    x_dat(ref_ind, 2) = -x_dat(ref_ind, 2)-2*br;
    
    ref_ind = x_dat(:, 3) > bl;
    x_dat(ref_ind, 3) = 2*bl - x_dat(ref_ind, 3);
    ref_ind = x_dat(:, 3) < -bl;
    x_dat(ref_ind, 3) = -x_dat(ref_ind, 3)-2*br;
    
    rho_dat = x_dat(:,1).^2/r2 + x_dat(:,2).^2/r2 + x_dat(:,3).^2/l2;
    intensity_accurate(j) = sum(exp(-2*rho_dat));
    intensity_efficient(j) = sum(exp(-2*rho_dat(inside_status)));
    
    if mod(j, m) == 1       
        % update inside or outside status
        inside_status = rho_dat < rho0;
    end % if mod==1
end % for j
end % function