clear;
clc;

% Iteration and Steps
h = 0.1;
N = 400/h;

% Initial Values
Vr = 2; % m3
na = zeros(1,N+1);
nb = zeros(1,N+1);
nc = zeros(1,N+1);
nd = zeros(1,N+1);
ne = zeros(1,N+1);
Tr = zeros(1,N+1);
t = zeros(1,N+1);

% Given Parameters
na(1) = 1000; nb(1) = 900; nc(1) = 150; nd(1) = 20; ne(1) = 0; Tr(1) = 330;
Kref = [5e-5, 2e-6];
Ea = [42500, 44500];
Dh = [-85000, -80000];
rho_m = 1100; % Kg/m3
Cpm = 4000; % J/(Kg K)
U = 5000; % J/(s m2 K)
A = 20; % m2
fluid_T = 340; % K
R = 8.314; % J/(mol K)

% Define derivatives
dnadt = @(r1, r2) -r1 - r2;
dnbdt = @(r1, r2) -r1;
dncdt = @(r1, r2) r1 - 2*r2;
dnddt = @(r1, r2) r2;
dnedt = @(r1, r2) r2;
dTdt = @(Qr, Qe) (Qr - Qe) / (rho_m * Cpm * Vr); % Added Vr here (total mass correction)

for n = 1:N
    t(n+1) = t(n) + h;
    
    % Step 1: Compute k1 (initial rates)
    K1 = Kref(1) * exp(-Ea(1)/R * (1/Tr(n) - 1/fluid_T));
    K2 = Kref(2) * exp(-Ea(2)/R * (1/Tr(n) - 1/fluid_T));
    r1 = K1 * (na(n)/Vr) * (nb(n)/Vr);
    r2 = K2 * (na(n)/Vr) * (nc(n)/Vr)^2;
    Qr = -(r1 * Dh(1) + r2 * Dh(2));
    Qe = U * A * (Tr(n) - fluid_T);
    
    k1a = h * dnadt(r1, r2);
    k1b = h * dnbdt(r1, r2);
    k1c = h * dncdt(r1, r2);
    k1d = h * dnddt(r1, r2);
    k1e = h * dnedt(r1, r2);
    k1T = h * dTdt(Qr, Qe);
    
    % Step 2: Compute k2 (midpoint)
    K1_mid = Kref(1) * exp(-Ea(1)/R * (1/(Tr(n)+k1T/2) - 1/fluid_T));
    K2_mid = Kref(2) * exp(-Ea(2)/R * (1/(Tr(n)+k1T/2) - 1/fluid_T));
    r1_mid = K1_mid * ((na(n)+k1a/2)/Vr) * ((nb(n)+k1b/2)/Vr);
    r2_mid = K2_mid * ((na(n)+k1a/2)/Vr) * ((nc(n)+k1c/2)/Vr)^2;
    Qr_mid = -(r1_mid * Dh(1) + r2_mid * Dh(2));
    Qe_mid = U * A * ((Tr(n)+k1T/2) - fluid_T);
    
    k2a = h * dnadt(r1_mid, r2_mid);
    k2b = h * dnbdt(r1_mid, r2_mid);
    k2c = h * dncdt(r1_mid, r2_mid);
    k2d = h * dnddt(r1_mid, r2_mid);
    k2e = h * dnedt(r1_mid, r2_mid);
    k2T = h * dTdt(Qr_mid, Qe_mid);
    
    % Step 3: Compute k3 (another midpoint)
    K1_mid2 = Kref(1) * exp(-Ea(1)/R * (1/(Tr(n)+k2T/2) - 1/fluid_T));
    K2_mid2 = Kref(2) * exp(-Ea(2)/R * (1/(Tr(n)+k2T/2) - 1/fluid_T));
    r1_mid2 = K1_mid2 * ((na(n)+k2a/2)/Vr) * ((nb(n)+k2b/2)/Vr);
    r2_mid2 = K2_mid2 * ((na(n)+k2a/2)/Vr) * ((nc(n)+k2c/2)/Vr)^2;
    Qr_mid2 = -(r1_mid2 * Dh(1) + r2_mid2 * Dh(2));
    Qe_mid2 = U * A * ((Tr(n)+k2T/2) - fluid_T);
    
    k3a = h * dnadt(r1_mid2, r2_mid2);
    k3b = h * dnbdt(r1_mid2, r2_mid2);
    k3c = h * dncdt(r1_mid2, r2_mid2);
    k3d = h * dnddt(r1_mid2, r2_mid2);
    k3e = h * dnedt(r1_mid2, r2_mid2);
    k3T = h * dTdt(Qr_mid2, Qe_mid2);
    
    % Step 4: Compute k4 (endpoint)
    K1_end = Kref(1) * exp(-Ea(1)/R * (1/(Tr(n)+k3T) - 1/fluid_T));
    K2_end = Kref(2) * exp(-Ea(2)/R * (1/(Tr(n)+k3T) - 1/fluid_T));
    r1_end = K1_end * ((na(n)+k3a)/Vr) * ((nb(n)+k3b)/Vr);
    r2_end = K2_end * ((na(n)+k3a)/Vr) * ((nc(n)+k3c)/Vr)^2;
    Qr_end = -(r1_end * Dh(1) + r2_end * Dh(2));
    Qe_end = U * A * ((Tr(n)+k3T) - fluid_T);
    
    k4a = h * dnadt(r1_end, r2_end);
    k4b = h * dnbdt(r1_end, r2_end);
    k4c = h * dncdt(r1_end, r2_end);
    k4d = h * dnddt(r1_end, r2_end);
    k4e = h * dnedt(r1_end, r2_end);
    k4T = h * dTdt(Qr_end, Qe_end);
    
    % Update variables
    na(n+1) = na(n) + (k1a + 2*k2a + 2*k3a + k4a)/6;
    nb(n+1) = nb(n) + (k1b + 2*k2b + 2*k3b + k4b)/6;
    nc(n+1) = nc(n) + (k1c + 2*k2c + 2*k3c + k4c)/6;
    nd(n+1) = nd(n) + (k1d + 2*k2d + 2*k3d + k4d)/6;
    ne(n+1) = ne(n) + (k1e + 2*k2e + 2*k3e + k4e)/6;
    Tr(n+1) = Tr(n) + (k1T + 2*k2T + 2*k3T + k4T)/6;
end

% Plotting
subplot(2,1,1);
plot(t, na, 'b', t, nb, 'r', t, nc, 'g', t, nd, 'm', t, ne, 'k');
xlabel('Time (s)');
ylabel('Moles (mol)');
legend('n_A', 'n_B', 'n_C', 'n_D', 'n_E');
grid on;

subplot(2,1,2);
plot(t, Tr, 'b');
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('T_r');
grid on;