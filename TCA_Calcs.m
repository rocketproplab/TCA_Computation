clear all; clc; format compact;
%% Thrust Chamber Design Parameters
% Using P_exit = 1, and plots from
% http://www.braeunig.us/space/comb-OM.htm P_Chamber, T_Chamber, y, Gas_Molecular_Weight can be extracted
%     after choosing the mixture ratio.
Characteristic_Length = 40;    % units - inches; interpolated use table 4-1 from Huzel+Huang Modern Engineering for Design of LRES 
Thrust = 2000;                 % units - pounds
T_Chamber_C = 2900.3;         % Chamber temperature; units - Celcius
P_Chamber = 20;                % Chamber pressure; units - atm
P_Exit = 0.82206329;                    % Exit pressure - sea level; units - atm
P_Ambient = .9278;                 % Ambient pressure - sea level; units - atm
Mixture_Ratio = 2.62;       % O/F mass ratio; unitless
Gas_Molecular_Weight = 18.833; % units - amu
y = 1.1510;                    % Ratio of Specific Heats - unitless
Working_Stress = 80000;        % units - psi; yield stress of inconel walling

%% Constants
g_earth =  32.1740;   % units - ft/seconds^2
R = 1545.348963;      % universal gas constant; units - ft*lbf/R*lb-mol
atm_to_psi = 14.6959; % units - psi/atm
pi = 3.141592654;
Factor_of_Safety = 1.4;

%% Adiabatic Flame Temperature - Chamber Temperature; CONVERSIONS
T_Chamber_K = T_Chamber_C+273.15; %units - Kelvin
T_Chamber_F = T_Chamber_C*1.8+32; %units - Fahrenheit
T_Chamber_R = T_Chamber_F+459.67; %units - Rankine

%% Expansion Ratio
Expansion_Ratio = (2/(y+1))^(1/(y-1)) * (P_Chamber/P_Exit)^(1/y)...
    / ( (y+1)/(y-1) * (1 - (P_Exit/P_Chamber)^(1-1/y)) )^.5 % unitless

%% Characteristic Velocity
c_star = (g_earth*y*R*T_Chamber_R/Gas_Molecular_Weight)^(.5)/(y*((2/(y+1))^((y+1)/(y-1)))^.5) % ft/second

%% Correction Factor
lamda = .985

%% Thrust Coefficient
Cf = lamda*( (2*y*y)/(y-1) * (2/(y+1))^((y+1)/(y-1)) * (1 - (P_Exit/P_Chamber)^((y-1)/y)) )^(1/2) + Expansion_Ratio*(P_Exit-P_Ambient)/P_Chamber

%% Specific Impulse
Isp = c_star*Cf/g_earth % units - seconds

%% Massflow
Massflow_Total = Thrust/Isp                        % units - lbs/s
Massflow_Fuel = Massflow_Total/(1+Mixture_Ratio)   % units - lbs/s
Massflow_Oxidizer = Massflow_Total - Massflow_Fuel % units - lbs/s

%% Nozzle Throat Parameters
T_Throat_R = T_Chamber_R*(2/(1+y))                            % units - Rankine
P_Throat = P_Chamber*(y/2 +.5)^(y/(1-y))                      % units - atm
Area_Throat = (Massflow_Total*Isp) / (P_Chamber*atm_to_psi*Cf)% units - in^2
Diameter_Throat = (Area_Throat/pi)^.5 * 2            % units - in

%% Contraction Ratio
Contraction_Ratio = 1.317303 + (32346880 - 1.317303)/(1 + (Diameter_Throat/1.041567e-10)^0.6774809) % units - inches;
% pulled from figure 4-9 of Modern Engineering for Design of LRES 
% fitted using https://mycurvefit.com/
%% Nozzle Exit
Area_Exit = Expansion_Ratio * Area_Throat % units - in^2
Diameter_Exit = (Area_Exit/pi)^.5 * 2 % units - in

%% Chamber Parameters
Volume_Chamber = Area_Throat*Characteristic_Length % units - in^3
Area_Chamber = Area_Throat*Contraction_Ratio       % units - in^2
Diameter_Chamber = (Area_Chamber/pi)^.5 * 2        % units - in^1
Length_Chamber = Volume_Chamber/Area_Chamber/1.1   % units - in^1

%% Wall Thickness
Thickness_ChamberWall = Factor_of_Safety*P_Chamber*atm_to_psi*Diameter_Chamber / (2*Working_Stress - 2*Factor_of_Safety*P_Chamber*atm_to_psi) %units - in^1

%% Nozzle Design
% y = Px+Q+(S*x+T)^(1/2)
% Curve Fitting Conditions  
%     1) y(x0 == 0) == y0 == 0
%     2) y(x1 == Nozzle_Length  - .382*Radius_Throat*sin(Theta_Initial) )
%     == y1 == Radius_exit - (.382*Radius_Throat*(1-cos(Theta_Initial)) +
%     Radius_Throat)
%     3) y'(x0 == 0) == tan(theta_initial)
%     4) y'(x1 == Nozzle_Length  - .382*Radius_Throat*sin(Theta_Initial)) == tan(theta_exit)
 
% http://soliton.ae.gatech.edu/people/jseitzma/classes/ae6450/bell_nozzle.pdf
% Sections 3.4+ in https://www.ewp.rpi.edu/hartford/~ernesto/S2013/EP/MaterialsforStudents/Lee/Sutton-Biblarz-Rocket_Propulsion_Elements.pdf
Radius_Throat = Diameter_Throat/2;
Radius_Exit = Diameter_Exit/2;

Theta_Initial = 25*3.141592/180; % Given in degrees by Rao, based on expansion ratio; here it's converted to radians
Theta_Exit = 12*3.141592/180; % Given in degrees by Rao, based on expansion ratio; here it's converted to radians

% (xN,yN) is the point where circular approximation ends and parabolic
% approximation begins; relative to center of nozzle throat
xN = .382*Radius_Throat*sin(Theta_Initial);                         %units inches
yN = .382*Radius_Throat*(1-cos(Theta_Initial)) + Radius_Throat;     %units inches

% (x1,y1) exit lip of the nozzle relative to (xN,yN)
x1 = .8*(Diameter_Throat/2)/tan(15*pi/180)*(Expansion_Ratio^.5-1+1.5*(1/cos(15*pi/180)-1)) - xN;%units inches
y1 = Radius_Exit -yN;                                           %units inches

N = tan(Theta_Initial); %Derivative at (x0,y0)
E = tan(Theta_Exit);    %Derivative at (x1,y1)

P = (E*(2*x1*N-y1) - y1*N)/(x1*N+x1*E-2*y1)
Q = -(y1-x1*E)^2*(x1*N-y1)/(2*(.5*x1*E+.5*x1*N-y1)^2)
S = -(y1-x1*E)^2*(y1-x1*N)^2*(E-N)/(2*(.5*x1*E+.5*x1*N-y1)^3)
T = (y1-x1*E)^4*(y1-x1*N)^2/(4*(.5*x1*E+.5*x1*N-y1)^4)
