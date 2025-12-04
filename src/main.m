%% Defining Constants and Equations
% Constants
params.hbar   = 1.0545718e-34;
params.m0     = 9.10938356e-31;
params.eV     = 1.60218e-19;
params.a      = 5.65325e-10;
params.Delta  = 0.34 * params.eV;
params.gamma1 = 6.8;
params.gamma2 = 2.1;
params.gamma3 = 2.73;

V_BZ = (2*pi/params.a)^3;
delta_E = 0.005 * params.eV;
E_min = -1.5 * params.eV;
E_max = 0;
E_bins = E_min:delta_E:E_max;
DOS_HH = zeros(size(E_bins));
DOS_LH = zeros(size(E_bins));

k_max = 0.025*2*pi/params.a;
resolution = 100;
k_range = linspace(-k_max , k_max, resolution);
E = zeros(6, length(k_range));

Q = @(kx,ky,kz) (params.hbar^2*params.gamma2/(2*params.m0)) * (kx^2+ky^2-2*kz^2);
R = @(kx,ky) (params.hbar^2/(2*params.m0)) * (-sqrt(3)*params.gamma2*(kx^2-ky^2) + 1i*2*sqrt(3)*params.gamma3*kx*ky);
S = @(kx,ky,kz) (params.hbar^2*params.gamma3/params.m0) * sqrt(3)*(kx-1i*ky)*kz;
P = @(kx,ky,kz) (params.hbar^2*params.gamma1/(2*params.m0)) * (kx^2+ky^2+kz^2);


%% Defining Matrix
% Defining Hamiltonian
for i = 1:length(k_range)
    % Define k
%{ 
    %For Problem 1
    kx = k_range(i); 
    ky = 0; 
    kz = 0; 
%}
    % For Problem 2
    kx = k_range(i)/sqrt(3); 
    ky = k_range(i)/sqrt(3); 
    kz = k_range(i)/sqrt(3);

    H = buildHamiltonian(kx, ky, kz, params);

    E(:,i) = eig(H);
end


%% Defining 3D Iso Surface Hamiltonian (Problem 3)
[kxx, kyy, kzz] = meshgrid(k_range, k_range, k_range);
E_HH = zeros(size(kxx));
E_LH = zeros(size(kxx));

for i = 1:length(k_range)
    for j = 1:length(k_range)
        for l = 1:length(k_range)
            kx_val = kxx(i,j,l)/sqrt(3);
            ky_val = kyy(i,j,l)/sqrt(3);
            kz_val = kzz(i,j,l)/sqrt(3);

            H = buildHamiltonian(kx_val, ky_val, kz_val, params);

            energies = eig(H);
            energies = sort(real(energies));

            E_HH(i,j,l) = energies(6);
            E_LH(i,j,l) = energies(4);
        end
    end
end

E_HH = E_HH / params.eV;
E_LH = E_LH / params.eV;
iso_energy_HH = -0.0004;
iso_energy_LH = -0.004;


%% Defining Density of States of Energy (DOS)
N_k = 1e5;

for i = 1:N_k
    total=2;
    while total>1.5
        randvar1 = (rand()*2)-1;
        randvar2 = (rand()*2)-1;
        randvar3 = (rand()*2)-1;

        total = abs(randvar1)+abs(randvar2)+abs(randvar3);
    end

    kx_4 = randvar1*2*pi/params.a;
    ky_4 = randvar2*2*pi/params.a;
    kz_4 = randvar3*2*pi/params.a;

    H = buildHamiltonian(kx_4, ky_4, kz_4, params);

    eigenvalues = sort(real(eig(H)));

    for j = 1:2
        if j==1
            band_Index = 6;
            E_n = eigenvalues(band_Index);
            if E_min <= E_n && E_n <= E_max
                bin_index = floor((E_n-E_min)/delta_E)+1;
                DOS_HH(bin_index) = DOS_HH(bin_index)+1;
            end
        elseif j==2
            band_Index = 4;
            E_n = eigenvalues(band_Index);
            if E_min <= E_n && E_n <= E_max
                bin_index = floor((E_n-E_min)/delta_E)+1;
                DOS_LH(bin_index) = DOS_LH(bin_index)+1;
            end
        end
    end
end

a_cm = params.a*100;
volume_factor = 4/(a_cm^3);
DOS_HH = DOS_HH/(N_k*delta_E/params.eV)*volume_factor;
DOS_LH = DOS_LH/(N_k*delta_E/params.eV)*volume_factor;
E_plot = E_bins / params.eV;


%% Trend Line for sqrt E
E_edge = min(E_plot);
DOS_trend = zeros(size(E_plot));
for i = 1:length(E_plot)
    if E_plot(i) >= E_edge
        DOS_trend(i) = sqrt(-1*E_plot(i));
    end
end
C = max(DOS_HH) / max(DOS_trend);
DOS_trend = C * DOS_trend;


%% Energy values for various wave vectors k (Problem 5)
k_points = [
    0.2, 0,   0;
    0.4, 0,   0;
    0.2, 0.2, 0
] * (2*pi/params.a);

k_labels = {'k = (0.2, 0, 0)', 'k = (0.4, 0, 0)', 'k = (0.2, 0.2, 0)'};

fprintf('Energy Results at Specified k-Points:\n');
fprintf('---------------------------------------\n');

for index = 1:size(k_points,1)
    kx = k_points(index,1);
    ky = k_points(index,2);
    kz = k_points(index,3);

    H = buildHamiltonian(kx, ky, kz, params);

    energies = sort(real(eig(H)));

    fprintf('%s:\n', k_labels{index});
    fprintf('   Heavy Hole (HH): %.6f eV\n', energies(6) / params.eV);
    fprintf('   Light Hole (LH): %.6f eV\n', energies(4) / params.eV);
    fprintf('   Split-off (SO):  %.6f eV\n', energies(2) / params.eV);
end


%% Deterministic DOS Calculation (Bonus)
DOS_deterministic = zeros(size(E_bins));
resolutionDOS = 50;
kx_range = linspace(-pi/params.a,pi/params.a,resolutionDOS);
ky_range = linspace(-pi/params.a,pi/params.a,resolutionDOS);
kz_range = linspace(-pi/params.a,pi/params.a,resolutionDOS);

for ix = 1:length(kx_range)
    for iy = 1:length(ky_range)
        for iz = 1:length(kz_range)
            kx = kx_range(ix);
            ky = ky_range(iy);
            kz = kz_range(iz);

            H = buildHamiltonian(kx, ky, kz, params);

            energies = sort(real(eig(H)));
            E_Deterministic = energies(6);

            if E_min <= E_Deterministic && E_Deterministic <= E_max
                bin_index = floor((E_Deterministic - E_min) / delta_E) + 1;
                DOS_deterministic(bin_index) = DOS_deterministic(bin_index) + 1;
            end
        end
    end
end

N_k_total = resolutionDOS^3;
DOS_deterministic = DOS_deterministic/(N_k_total*delta_E/params.eV)/(a_cm^3);


%% Graphing Results
energies_eV = E / params.eV;

figure;
hold on;
for j=1:6
    plot(k_range, energies_eV(j,:),'LineWidth',1.5);
end
hold off;
xlabel('Wavevector k (1/m)');
ylabel('Energy (eV)');
title('Electronic Band Structure (GaAs)');
legend('Band 1','Band 2','Band 3','Band 4','Band 5','Band 6');
grid on;

figure;
isosurface(kxx, kyy, kzz, E_HH, iso_energy_HH);
FaceColor = 'red';
EdgeColor = 'none';
title('3D Iso-Energy Surface: Heavy Hole');
xlabel('k_x (1/m)'); ylabel('k_y (1/m)'); zlabel('k_z (1/m)');
grid on;
axis equal;
camlight;
lighting gouraud;

figure;
isosurface(kxx, kyy, kzz, E_LH, iso_energy_LH);
title('3D Iso-Energy Surface: Light Hole');
xlabel('k_x (1/m)'); ylabel('k_y (1/m)'); zlabel('k_z (1/m)');
grid on;
axis equal;
camlight;
lighting gouraud;

figure;
plot(E_plot, DOS_HH, E_plot, DOS_LH, E_plot, DOS_trend, 'LineWidth', 1.5);
hold on;
plot(E_bins / params.eV, DOS_deterministic, 'LineWidth', 2);
xlabel('Energy (eV)');
ylabel('Density of States (eV^{-1} cm^{-3})');
title('Density of States for Heavy Hole and Light Hole Bands');
legend('Heavy hole','Light Hole','Square Root of Energy Trend', 'Deterministic DOS');
grid on;
hold off;