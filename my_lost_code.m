clear;
clc;
close all;

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
R = 8314; % (J/mol K) 

dTdt = @(~, Qr, Qe) (Qr - Qe)/(rho_m * Cpm);
dnadt = @(r1, r2) Vr - r1 - r2;
dnbdt = @(r1, r2) Vr - r1;
dncdt = @(r1, r2) Vr + r1 - 2*r2;
dnddt = @(r1, r2) Vr + r2;
dnedt = @(r1, r2) Vr + r2;


for n = 1:N
    t(n+1) = t(n) + h;

    K1(n) = Kref(1) * exp((Ea(1)/R) * ((1/323)-(1/Tr(n))));
    K2(n) = Kref(2) * exp((Ea(2)/R) * ((1/323)-(1/Tr(n))));

    r1(n) = K1(n) * (na(n)/Vr) * (nb(n)/Vr);
    r2(n) = K2(n) * (na(n)/Vr) * (nc(n)/Vr)^2;

    Qe(n) = U * A * (Tr(n) - fluid_T)/Vr;
    Qr(n) = -(r1(n) * Dh(1)) + (r2(n) * Dh(2));

    K1T_ = h * dTdt(t(n),Qr(n), Qe(n));
    k1a_ = h * dnadt(r1(n),r2(n));
    k1b_ = h * dnbdt(r1(n), r2(n));
    k1c_ = h * dncdt(r1(n), r2(n));
    k1d_ = h * dnddt(r1(n), r2(n));
    k1e_ = h * dnedt(r1(n), r2(n));

    K2T_ = h * dTdt(t(n) + h/2, Qr(n) + h/2, Qe(n) + h/2);
    k2a_ = h * dnadt(r1(n) + k1a_/2, r2(n) + k1b_/2);
    k2b_ = h * dnbdt(r1(n) + k1a_/2, r2(n) + k1b_/2);
    k2c_ = h * dncdt(r1(n) + k1a_/2, r2(n) + k1b_/2);
    k2d_ = h * dnddt(r1(n) + k1a_/2, r2(n) + k1b_/2);
    k2e_ = h * dnedt(r1(n) + k1a_/2, r2(n) + k1b_/2);

    K3T_ = h * dTdt(t(n) + h/2, Qr(n) + h/2, Qe(n) + h/2);
    k3a_ = h * dnadt(r1(n) + k2a_/2, r2(n) + k2b_/2);
    k3b_ = h * dnadt(r1(n) + k2a_/2, r2(n) + k2b_/2);
    k3c_ = h * dnadt(r1(n) + k2a_/2, r2(n) + k2b_/2);
    k3d_ = h * dnadt(r1(n) + k2a_/2, r2(n) + k2b_/2);
    k3e_ = h * dnadt(r1(n) + k2a_/2, r2(n) + k2b_/2);

    K4T_ = h * dTdt(t(n) + h/2, Qr(n) + h, Qe(n) + h);
    k4a_ = h * dnadt(r1(n) + k3a_, r2(n) + k3b_);
    k4b_ = h * dnadt(r1(n) + k3a_, r2(n) + k3b_);
    k4c_ = h * dnadt(r1(n) + k3a_, r2(n) + k3b_);
    k4d_ = h * dnadt(r1(n) + k3a_, r2(n) + k3b_);
    k4e_ = h * dnadt(r1(n) + k3a_, r2(n) + k3b_);

    %Update Tr
    Tr(n+1) = Tr(n) + ((1/6)*( K1T_ + 2*K2T_ + 2*K3T_ + K4T_ ));
    na(n+1) = na(n) + ((1/6)*( k1a_ + 2*k2a_ + 2*k3a_ + k4a_ ));
    nb(n+1) = nb(n) + ((1/6)*( k1b_ + 2*k2b_ + 2*k3b_ + k4b_ ));
    nc(n+1) = nc(n) + ((1/6)*( k1c_ + 2*k2c_ + 2*k3c_ + k4c_ ));
    nd(n+1) = nd(n) + ((1/6)*( k1d_ + 2*k2d_ + 2*k3d_ + k4d_ ));
    ne(n+1) = ne(n) + ((1/6)*( k1e_ + 2*k2e_ + 2*k3e_ + k4e_ ));

end
 plot(t, Tr, 'b');
 grid on;
