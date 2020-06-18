%% Simulate and Fit for a batch of trajectories
%% Simulate and analyze one trajectory at a time,
%% write the data into a files,
%% perform least squares fit for five methods,
%% write the fitted data to a file.

tic;

%% Write the parameters to an info file
filename = 'sim-info.txt';
fileID = fopen(filename,'w');

ns = 200;
fprintf(fileID, "-Number of simulations: %d.\n", ns);

dt = 1e-06;
fprintf(fileID,"Time step dt: %g.\n", dt);


T_dat = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01];
nt_dat = round(T_dat/dt);
%nt_dat = int64(nt_dat);

nT = length(T_dat);
fprintf(fileID, "T for temporal average are: \n");
for i = 1:nT
    fprintf(fileID, "%g, ", T_dat(i));
end
fprintf(fileID, "\n");
Tmax = max(T_dat);


taumax = 2e-04;
fprintf(fileID, "Tau max for correlation: %g.\n", taumax);
ntau = floor(taumax/dt);
%ntau = int64(ntau);

D = 1.3e-10;
fprintf(fileID, "-Diffusivity: %g.\n", D);

r = 2e-7;
fprintf(fileID, "Optical system r: %g.\n", r);
l = 1e-6;
fprintf(fileID, "Optical system l: %g.\n", l);

br = 8*r;
fprintf(fileID, "Box length: %g.\n", 2*br);
bl = 8*l;
fprintf(fileID,"Box height: %g.\n", 2*bl);
conc = 0.02;
fprintf(fileID, "Concentration: %g.\n", conc);

% derived simulation parameters 
vol_sol = 8*br*br*bl;
avogadro = 6.0221409e+23;
npart_d = conc*avogadro*(1e-06)*1000*vol_sol;
npart_i = round(npart_d);
%npart_i = int64(npart_i);
fprintf(fileID, "Number of particles %d.\n", npart_i);

vol_det = pi * sqrt(pi) * r*r*l;
neff = conc*avogadro*(1e-06)*1000*vol_det;
fprintf(fileID, "-Effective number of particle in illuminatd region: %g.\n", neff);

% ntime: the number of time points in a trajectory
ntime = round((Tmax+taumax)/dt); 
%fprintf(fileID, "-The number of time points: %d.\n", ntime);
%nt = ntime - ntau;
t_dat = (1:ntime)'*dt;
tau_dat = (1:ntau)'*dt;

% efficient simulation
a = 3*r;
c = 3*l;
m = 20;
fprintf(fileID, "Efficient simulation important region: (x1/a)^2 + (x2/a)^2 + (x3/c)^2 <= 1.\n");
fprintf(fileID, "a = %g.\n", a);
fprintf(fileID, "c = %g.\n", c);
fprintf(fileID, "m = %g.\n", m);

fclose(fileID); % close the info file

%% Simulate the first trajectory and write into binary files

% simulate intensity trajectory, create a file and write
intensity = get_intensity_eff(dt, ntime, D, br, bl, npart_i, r, l, m);

write_intensity_status = 0;
write_corr_status = 0;

if write_intensity_status == 1
intensity_filename = 'intensity1.bin';
fileID = fopen(intensity_filename,'w');
fwrite(fileID,intensity,'double');
fclose(fileID);
end


%%
% Corr, Corr3, G, g
Corr = zeros(ntau, 1);
Corr3 = zeros(ntau, 1);
G = zeros(ntau, 1);
g2 = zeros(ntau, 1);

corr_filenames = ["corr_T1.bin", "corr_Tp5.bin", "corr_Tp2.bin",...
    "corr_Tp1.bin", "corr_Tp05.bin", "corr_Tp02.bin", "corr_Tp01.bin"];
corr3_filenames = ["corr3_T1.bin", "corr3_Tp5.bin", "corr3_Tp2.bin",...
    "corr3_Tp1.bin", "corr3_Tp05.bin", "corr3_Tp02.bin", "corr3_Tp01.bin"];
G_filenames = ["G_T1.bin", "G_Tp5.bin", "G_Tp2.bin",...
    "G_Tp1.bin", "G_Tp05.bin", "G_Tp02.bin", "G_Tp01.bin"];
g2_filenames = ["g2_T1.bin", "g2_Tp5.bin", "g2_Tp2.bin",...
    "g2_Tp1.bin", "g2_Tp05.bin", "g2_Tp02.bin", "g2_Tp01.bin"];

D_corr_filenames = ["D_corr_T1.bin", "D_corr_Tp5.bin", "D_corr_Tp2.bin",...
    "D_corr_Tp1.bin", "D_corr_Tp05.bin", "D_corr_Tp02.bin", "D_corr_Tp01.bin"];
D_corr3_filenames = ["D_corr3_T1.bin", "D_corr3_Tp5.bin", "D_corr3_Tp2.bin",...
    "D_corr3_Tp1.bin", "D_corr3_Tp05.bin", "D_corr3_Tp02.bin", "D_corr3_Tp01.bin"];
D_G_filenames = ["D_G_T1.bin", "D_G_Tp5.bin", "D_G_Tp2.bin",...
    "D_G_Tp1.bin", "D_G_Tp05.bin", "D_G_Tp02.bin", "D_G_Tp01.bin"];
D_g2_filenames = ["D_g2_T1.bin", "D_g2_Tp5.bin", "D_g2_Tp2.bin",...
    "D_g2_Tp1.bin", "D_g2_Tp05.bin", "D_g2_Tp02.bin", "D_g2_Tp01.bin"];
D_2n3_filenames = ["D_2n3_T1.bin", "D_2n3_Tp5.bin", "D_2n3_Tp2.bin",...
    "D_2n3_Tp1.bin", "D_2n3_Tp05.bin", "D_2n3_Tp02.bin", "D_2n3_Tp01.bin"];

Ne_ig = zeros(nT, 1); 
F0_ig = zeros(nT, 1); 


% Models 
kd = 10^10;
% Gfun_model is the scaled autocorrelation function 
% since D is on the 10^(-10) scale while Neff is on the scale of 1
Gfun_model = @(params, tau_dat) 1/params(1) ...
    * (1./(1+4*params(2)*tau_dat/(kd*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(kd*l^2)));
% where params(1) is Neff, the effective number of particles in illuminated volume
% and params(2) is D, the diffusivity

Corr_model = @(params, tau_dat) params(3)^2 * ( params(1) * ...
    (1/8./(1+4*params(2)*tau_dat/(kd*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(kd*l^2))) ...
    + params(1)^2/8);
% where params(1) is Neff, parmas(2) is D, and params(3) is F0

Corr3_model = @(params, tau_dat) params(3)^3 * ( params(1)/sqrt(8) * ...
    1./(3+16*params(2)*tau_dat/(kd*r^2)) .* 1./sqrt(3+16*params(2)*tau_dat/(kd*l^2)) ...
    + params(1)^3/(8*sqrt(8)) + params(1)^2/(8*sqrt(8))...
    + params(1)^2/(8*sqrt(8))*2./(1+4*params(2)*tau_dat/(kd*r^2)) ...
    .* (1./sqrt(1+4*params(2)*tau_dat/(kd*l^2))) );
% where params(1) is Neff, parmas(2) is D, and params(3) is F0

Corr_2n3_model = @(params, tau_dat) ...
    Corr_model(params, tau_dat) + Corr3_model(params, tau_dat);

M2bar_model = @(param, tau_dat) ( param(1) * ...
    (1/8./(1+4*param(2)*tau_dat/(kd*r^2))) .* (1./sqrt(1+4*param(2)*tau_dat/(kd*l^2))) ...
    + param(1)^2/8);
M3bar_model = @(param, tau_dat) ( param(1)/sqrt(8) * ...
    1./(3+16*param(2)*tau_dat/(kd*r^2)) .* 1./sqrt(3+16*param(2)*tau_dat/(kd*l^2)) ...
    + param(1)^3/(8*sqrt(8)) + param(1)^2/(8*sqrt(8))...
    + param(1)^2/(8*sqrt(8))*2./(1+4*param(2)*tau_dat/(kd*r^2)) ...
    .* (1./sqrt(1+4*param(2)*tau_dat/(kd*l^2))) );

g2_model = @(params, tau_dat) (1./(1+4*params*tau_dat/(kd*r^2))) ...
    .* (1./sqrt(1+4*params*tau_dat/(kd*l^2)));

%%


for i=1:nT
    nt = nt_dat(i);
    I_1 = mean(intensity(1:nt));
    for k = 1:ntau
        Corr(k) = intensity(1:nt)'*intensity(k:(k+nt-1))/nt;
        Corr3(k) = intensity(1:nt).^2'*intensity(k:(k+nt-1))/nt;
        I_k = mean(intensity(k: (k+nt-1)));
        G(k) = Corr(k)/(I_1*I_k) - 1;
        g2(k) =  (Corr(k) - I_1*I_k)/(std(intensity(1:nt))*std(intensity(k:(k+nt-1))));
    end %for tau_k
    
    Ne_ig(i) = 1/G(1);
    F0_ig(i) = mean(intensity(1:nt))*sqrt(8)*G(1);
    
    
    % -----------fit Corr
    params0 = [Ne_ig(i), 1, F0_ig(i)];
    % f1 maps params to ntau dimensional vector
    f1 = @(params)  Corr - Corr_model(params, tau_dat);
    options = optimoptions('lsqnonlin','Display', 'off', 'TolX', 1e-6);
    params_corr = lsqnonlin(f1, params0, [0,0,0], [], options);
    
    fileID = fopen(D_corr_filenames(i), 'w');
    fwrite(fileID, params_corr(2)/kd, 'double');
    fclose(fileID);
    
    % ------------fit Corr3
    params0 = [Ne_ig(i), 1, F0_ig(i)];
    f1 = @(params)  Corr3 - Corr3_model(params, tau_dat);
    params_corr3 = lsqnonlin(f1, params0, [0,0,0], [], options);
    
    fileID = fopen(D_corr3_filenames(i), 'w');
    fwrite(fileID, params_corr3(2)/kd, 'double');
    fclose(fileID);

    % ------------fit G
    params0 = [Ne_ig(i), 1];
    f1 = @(params)  G - Gfun_model(params, tau_dat);
    params_G = lsqnonlin(f1, params0, [0, 0], [], options);
    
    fileID = fopen(D_G_filenames(i), 'w');
    fwrite(fileID, params_G(2)/kd, 'double');
    fclose(fileID);
    
    % ------------fit g2
    % f1 maps params to ntau dimensional vector
    f1 = @(params)  g2 - g2_model(params, tau_dat);
    D_g2 = lsqnonlin(f1, 1, 0, [], options);
    
    fileID = fopen(D_g2_filenames(i), 'w');
    fwrite(fileID, D_g2/kd, 'double');
    fclose(fileID);
    
    %-------------fit Corr and Corr3
    f2 = @(params) [Corr/params(3)^2 ...
        - M2bar_model([params(1), params(2)], tau_dat); ...
        Corr3/params(3)^3 ...
        - M3bar_model([params(1), params(2)], tau_dat)];
    params_2n3 = lsqnonlin(f2, params_corr, [0,0,0], [], options);

    fileID = fopen(D_2n3_filenames(i), 'w');
    fwrite(fileID, params_2n3(2)/kd, 'double');
    fclose(fileID);
    
    % create and write into data files
    if write_corr_status == 1
    fileID = fopen(corr_filenames(i),'w');
    fwrite(fileID, Corr, 'double');
    fclose(fileID);
    fileID = fopen(corr3_filenames(i),'w');
    fwrite(fileID, Corr3, 'double');
    fclose(fileID);
    fileID = fopen(G_filenames(i),'w');
    fwrite(fileID, G, 'double');
    fclose(fileID);
    fileID = fopen(g2_filenames(i),'w');
    fwrite(fileID, g2, 'double');
    fclose(fileID);
    end
end % for T_dat(i)


%% Simulation more trajectories and append the binary files

for j = 2:ns
    intensity = get_intensity_eff(dt, ntime, D, br, bl, npart_i, r, l, m);
    for i=1:nT
        nt = nt_dat(i);
        I_1 = mean(intensity(1:nt));
        for k = 1:ntau
            Corr(k) = intensity(1:nt)'*intensity(k:(k+nt-1))/nt;
            Corr3(k) = intensity(1:nt).^2'*intensity(k:(k+nt-1))/nt;
            I_k = mean(intensity(k: (k+nt-1)));
            G(k) = Corr(k)/(I_1*I_k) - 1;
            g2(k) =  (Corr(k) - I_1*I_k)/(std(intensity(1:nt))*std(intensity(k:(k+nt-1))));
        end %for tau_k
        
        Ne_ig(i) = 1/G(1);
        F0_ig(i) = mean(intensity(1:nt))*sqrt(8)*G(1);
        
        
        % -----------fit Corr
        params0 = [Ne_ig(i), 1, F0_ig(i)];
        % f1 maps params to ntau dimensional vector
        f1 = @(params)  Corr - Corr_model(params, tau_dat);
        options = optimoptions('lsqnonlin','Display', 'off', 'TolX', 1e-6);
        params_corr = lsqnonlin(f1, params0, [0,0,0], [], options);
        
        fileID = fopen(D_corr_filenames(i), 'a');
        fwrite(fileID, params_corr(2)/kd, 'double');
        fclose(fileID);
        
        % ------------fit Corr3
        params0 = [Ne_ig(i), 1, F0_ig(i)];
        f1 = @(params)  Corr3 - Corr3_model(params, tau_dat);
        params_corr3 = lsqnonlin(f1, params0, [0,0,0], [], options);
        
        fileID = fopen(D_corr3_filenames(i), 'a');
        fwrite(fileID, params_corr3(2)/kd, 'double');
        fclose(fileID);
        
        % ------------fit G
        params0 = [Ne_ig(i), 1];
        f1 = @(params)  G - Gfun_model(params, tau_dat);
        params_G = lsqnonlin(f1, params0, [0, 0], [], options);
        
        fileID = fopen(D_G_filenames(i), 'a');
        fwrite(fileID, params_G(2)/kd, 'double');
        fclose(fileID);
        
        % ------------fit g2
        % f1 maps params to ntau dimensional vector
        f1 = @(params)  g2 - g2_model(params, tau_dat);
        D_g2 = lsqnonlin(f1, 1, 0, [], options);
        
        fileID = fopen(D_g2_filenames(i), 'a');
        fwrite(fileID, D_g2/kd, 'double');
        fclose(fileID);
        
        %-------------fit Corr and Corr3
        f2 = @(params) [Corr/params(3)^2 ...
            - M2bar_model([params(1), params(2)], tau_dat); ...
            Corr3/params(3)^3 ...
            - M3bar_model([params(1), params(2)], tau_dat)];
        params_2n3 = lsqnonlin(f2, params_corr, [0,0,0], [], options);
        
        fileID = fopen(D_2n3_filenames(i), 'a');
        fwrite(fileID, params_2n3(2)/kd, 'double');
        fclose(fileID);
        
        % create and write into data files
        if write_corr_status == 1
            fileID = fopen(corr_filenames(i),'a');
            fwrite(fileID, Corr, 'double');
            fclose(fileID);
            fileID = fopen(corr3_filenames(i),'a');
            fwrite(fileID, Corr3, 'double');
            fclose(fileID);
            fileID = fopen(G_filenames(i),'a');
            fwrite(fileID, G, 'double');
            fclose(fileID);
            fileID = fopen(g2_filenames(i),'a');
            fwrite(fileID, g2, 'double');
            fclose(fileID);
        end
    end % for T_dat(i)
end

toc;