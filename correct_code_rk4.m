%% Question 2.6  Riccardo Tesser, Vincenzo Russo Advanced Reactor Modeling with MATLAB
%% Continuation in rxn_2.m
close all;
global Vr delta_H Ea U A rho H_av Ts kref Tref R Cpm
% Constants
R = 8.314; % J/mol.K
delta_H = [-85000 -80000]; % J/mol
Ea = [42500 44500];       % J/mol
% Initial Conditions
y = [1000 900 150 20 0 330]; % [A B C D E T]
Vr = 2;           % Reactor volume in L
U = 5000;         % Heat transfer coefficient, W/m^2.K
A = 20;           % Area for heat exchange, m^2
rho = 1100;       % kg/m^3
Cpm = 4000;        %  J/kg.K
Ts = 340;         % Surrounding temperature, K
kref = 1e-6 * [50 2]; % Reference rate constants
Tref = 323;       % Reference temperature, K
% Time span
h = 0.1;  
tspan = 1:h:400; 
N = length(tspan);  
% Solve ODE
%[t, y] = ode45(@reactorODE, tspan, mol_init);
for i = 1:N-1
    t = tspan(i);
    n = y(i,:)';
    k1 = h * reactorODE(t, n);
    k2 = h * reactorODE(t + h/2, n + k1/2);
    k3 = h * reactorODE(t + h/2, n + k2/2);
    k4 = h * reactorODE(t + h, n + k3);
    y(i+1,:) = (n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))';
end
figure;
subplot(2,1,1);  
plot(tspan, y(:,1:5), 'LineWidth', 1.5);
legend('A', 'B', 'C', 'D', 'E', 'Location', 'best');
xlabel('Time');
ylabel('Moles');
title('Mole Profiles in Reactor');
grid on;
subplot(2,1,2);  
plot(tspan, y(:,6), 'r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Temperature (K)');
title('Temperature Profile in Reactor');
grid on;
% -----------------------------
% ODE Function
function dydt = reactorODE(t, y)
    global Vr delta_H Ea U A rho Cpm Ts kref Tref R
  
    nA = y(1); nB = y(2); nC = y(3);
    nD = y(4); nE = y(5); T = y(6);
    % Concentrations (mol/L)
    Ca = nA / Vr;
    Cb = nB / Vr;
    Cc = nC / Vr;
    %Cd = nD / Vr;
    %Ce = nE / Vr;
 
    k = kref .* exp( Ea./R .* (1./Tref - 1./T) );
    % Reaction rates
    r1 = k(1) * Ca * Cb;
    r2 = k(2) * Ca * Cc^2;
    % Molar rate equations
    dnAdt = Vr * (-r1 - r2);
    dnBdt = Vr * (-r1);
    dnCdt = Vr * ( r1 - 2*r2 );
    dnDdt = Vr * ( r2 );
    dnEdt = Vr * ( r2 );
    % Energy balance
    Qr = -(r1 * delta_H(1) + r2 * delta_H(2)); % J/s
    QE = U * A * (T - Ts)/Vr;                          % J/s
    dTdt = (Qr - QE) / (rho * Cpm );           % K/s
    dydt = [dnAdt; dnBdt; dnCdt; dnDdt; dnEdt; dTdt];
end