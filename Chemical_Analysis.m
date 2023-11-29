% Matthew Simpson
% Propellant Analysis

%% GOX - Methane

clc;clear

mixture_ratio = 4;

GOX_MW = 2*16;
Methane_MW = 16;

a = 1; % number of moles of RP1
b = (mixture_ratio*Methane_MW)/GOX_MW;

LHS = [1 0 1; 0 2 4;2 1 0];
RHS = [1;4;b*2];

Sol = LHS\RHS;
c = Sol(1);
d = Sol(2);
e = Sol(3);

MW_CO2 = 44.009;
MW_H20 = 18.0146;
MW_Methane = 16.04;

Mass_CO2 = MW_CO2*c;
Mass_H20 = MW_H20*d;
Mass_Methane = MW_Methane*e;

MW_Mixture = (Mass_CO2+Mass_H20+Mass_Methane) / (c+d+e);

total_moles = c+d+e;
Mole_frac_CO2 = c/total_moles;
Mole_frac_H20 = d/total_moles;
Mole_frac_Methane = e/total_moles;

CP_H20 = 56.583;
CP_CO2 = 62.573;
CP_Methane = 35.706; %using the book table

RU = 8314.5;
Cp_Mixture = (CP_H20*Mole_frac_H20) + (CP_CO2*Mole_frac_CO2) + (CP_Methane*Mole_frac_Methane);
gamma_mix = Cp_Mixture / (Cp_Mixture -(RU/1000));

R = RU/MW_Mixture;


% compute the adibatic flame temperature assuming complete combustion since
% mixture ratio is stoichiometric

H0F_Methane = -17.8890*1000;
H0F_O2 = 0;
H0F_CO2 = -94.0518*1000;
H0F_H20 = -57.7979*1000;

syms Tc T
CP_CO2 = 75.513 - (0.18732*10^-3)*T - 661.85*T^(-.5);
CP_Water = 29.182 + (14.503*(T/1000)) + (-2.0235*(T/1000)^2);
equ1 = 1*(H0F_Methane)+b*(H0F_O2) == c*(H0F_CO2+int(CP_CO2,[298,Tc])) + d*(H0F_H20+int(CP_Water,[298,Tc]));

Adiabatic_flame_temp = vpasolve(equ1,Tc);

C_star = sqrt((R*Adiabatic_flame_temp)/gamma_mix)*((gamma_mix+1)/2)^((gamma_mix+1)/(2*(gamma_mix-1)));

% example values for the rocket (modify these)

A_throat = 0.01; %m^2
mdot = 5;
pe = 20000;
pa = 101325;
episolon = 5.6;

Pc = (C_star*mdot)/A_throat;

pt1 = (2*gamma_mix^2)/(gamma_mix-1);
pt2 = (2/(gamma_mix+1))^((gamma_mix+1)/(gamma_mix-1));
pt3 = (1-(pe/Pc))^((gamma_mix-1)/gamma_mix);

cf = (pt1*pt2*pt3)^(.5)+((pe/Pc)-(pa/Pc))*episolon;

ISP = cf*C_star/9.81;

fprintf('Mole Weight of Mixture: %.4f g/mol\n',MW_Mixture)
fprintf('CP of Mixture: %0.4f J/mol*k\n',Cp_Mixture)
fprintf('γ of Mixture: %0.4f\n',gamma_mix)
fprintf('Adiabatic Flame Temperature: %0.4f K\n',Adiabatic_flame_temp)
fprintf('C* : %0.4f m/s\n',C_star)
fprintf('Thrust Coefficient: %0.4f \n',cf)
fprintf('ISP: %0.4f \n',ISP)