s = struct();

s.r = 2;           % Reactor volume in L
s.U = 5000;         % Heat transfer coefficient, W/m^2.K
s.A = 20;           % Area for heat exchange, m^2
s.rho = 1100;       % kg/m^3
s.Cpm = 4000;        %  J/kg.K
s.Ts = 340;         % Surrounding temperature, K
s.kref = 1e-6 * [50 2]; % Reference rate constants
s.Tref = 323;       % Reference temperature, K

m =10
display(ode(7,4,6))
ode(1,4,m)


function param = ode(a,k,z)
    dkdt = z.r * 5 + 6;
   
end