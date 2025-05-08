clear;
clc;
% close all;

%Iteration and Steps
h = 0.1;
N = 400/h;


%Initial Values
Vr = 2; % M3
na = zeros(1,N);
nb = zeros(1,N);
nc = zeros(1,N);
nd = zeros(1,N);
ne = zeros(1,N);
Tr = zeros(1,N);
K1 = zeros(1,N);
K2 = zeros(1,N);
t = zeros(1,N);
r1 = zeros(1,N);
r2 = zeros(1,N);
Qe = zeros(1,N);
Qr = zeros(1,N);

%Given Parameters
na(1)= 1000; nb(1) = 900; nc(1) = 150; nd(1) = 20; ne(1) = 0; Tr(1) = 330;
Kref = [5e-5 2e-6];
Ea = [42500 44500];
Dh = [-85000 -80000];
rho_m = 1100; %Kg/m3
Cpm = 4000; %(J/(Kg K))
U = 5000; %(J/(s m2 K))
A = 20; %m2
fluid_T = 340; %K
R = 8.314; % (J/mol K) 

dTdt = @(Qr, Qe) (Qr - Qe)/(rho_m * Cpm);
dnadt = @(t, r1, r2) Vr * (- r1 - r2);
dnbdt = @(t, r1, r2) Vr * (-r1) ;
dncdt = @(t, r1, r2) Vr * (r1 - 2*r2);
dnddt = @(t, r1, r2) Vr * r2;
dnedt = @(t, r1, r2) Vr * r2;


for n = 1:N
    t(n+1) = t(n) + h;

    Ca(n) = na(n)/Vr;
    Cb(n) = nb(n)/Vr; 
    Cc(n) = nc(n)/Vr;

    k1(n) = Kref(1) * exp((Ea(1)/R) * ((1/323)-(1/Tr(n))));
    k2(n) = Kref(2) * exp((Ea(2)/R) * ((1/323)-(1/Tr(n))));

    r1(n) = k1(n) * Ca(n) * Cb(n);
    r2(n) = k2(n) * Ca(n) * Cc(n)^2;

    Qr(n) = -(r1(n) * Dh(1)) + (r2(n) * Dh(2));
    Qe(n) = U * A * (Tr(n) - fluid_T);

    k1_T = h * dTdt( Qr(n), Qe(n));
    k1_na = h * dnadt(t(n), r1(n), r2(n));
    k1_nb = h * dnbdt(t(n), r1(n), r2(n));
    k1_nc = h * dncdt(t(n), r1(n), r2(n));
    k1_nd = h * dnddt(t(n), r1(n), r2(n));
    k1_ne = h * dnedt(t(n), r1(n), r2(n));

    k2_T = h * dTdt(Qr(n) + h/2, Qe(n) + h/2);
    k2_na = h * dnadt(t(n) + k1_na/2, r1(n) + h/2, r2(n) + h/2);
    k2_nb = h * dnbdt(t(n) + k1_nb/2, r1(n) + h/2, r2(n) + h/2);
    k2_nc = h * dnbdt(t(n) + k1_nc/2, r1(n) + h/2, r2(n) + h/2);
    k2_nd = h * dnbdt(t(n) + k1_nd/2, r1(n) + h/2, r2(n) + h/2);
    k2_ne = h * dnbdt(t(n) + k1_ne/2, r1(n) + h/2, r2(n) + h/2);

    k3_T = h * dTdt( Qr(n) + h/2, Qe(n) + h/2);
    % k3_na = h * dnadt(t(n) + k2_na, r1(n) + h/2, )


    k4_T = h * dTdt( Qr(n) + h, Qe(n) + h);

    Tr(n+1) = Tr(n) + (1/6) * (k1_T + 2*k2_T + 2*k3_T + k4_T);

  
end
   plot(t, Tr, 'b')
% subplot(2,1,1); plot(t, na);
% 
% hold on;
% plot(t, nb);
% plot(t, nc);
% plot(t, nd);
% plot(t, ne);
% xlabel('Time');
% ylabel('Moles');
% % xlim([0 800])
% % ylim([0 1200])
% legend('n_A', 'n_B', 'n_C', 'n_D', 'n_E');
% grid on;
% 
% subplot(2,1,2); plot(t, Tr, 'b');
% legend('n_T');
% grid on;
% 
% hold off;