%% 102B Project SPRING 2021
%Name: Linh Tang



clc; clear all; close all;

global P A B C;


%% Question 2a
% a mixture composition of 17% toluene, 23% benzene, 29% cyclohexane, and
% 31% chlorobenzene

% Let this mixture be fully liquid. At a pressure of 1 atm, what is the bubble point
% temperature and composition of the first bubble?

% Assume ideal behavior

% P = 1atm = 101325 kPa
P = 101.325;
epsilon = 1e-6;
x_new = [0.17, 0.23, 0.29, 0.31];
n = length(x_new);

% Initialize 
phi (1:n) = 1;
gamma(1:n) = 1;


y_new = zeros(1,n);


T_sat = Cal_T (P,n);

% Calculate T = sum( xi*Ti_sat)(Eq 14.12)
T1 = x_new*T_sat';

%fprintf("T (eq.14.12: %d", T1);

% Calculate Psat
P_sat = Cal_P(T1, n);


% Choose toluene as species j
A= 13.9320; 
B= 3056.96; 
C= 217.625;

Pj_sat = P/((x_new.*gamma./phi./P_sat(1))*P_sat');

%Calculate T with A, B, C (Eq 14.15)
T2 = (B /(A - log(Pj_sat)))- C;

delta_T = abs(T1 - T2);


while (delta_T > epsilon)
    T1 = T2;
    P_sat = Cal_P(T1,n);
    
    y_new=x_new.*gamma.*P_sat./(P*phi);
    
    % Calculate new Pj_sat
    Pj_sat_new = exp(A-(B/(T1+C))); % from Antoine
    Pj_sat = P/((x_new./Pj_sat_new)*P_sat');
    
    %Calculate new T2
    T2 = (B/(A-log(Pj_sat)))-C; % using Antoine
    
    %calc new delta
    delta_T = abs(T2 - T1);
    
    
    
end

% Store values into variables
SOL_q2a_T = T2;
SOL_q2a_y = y_new;



%% Question 2b

% Assume Ideal behavior, fully gas, calculate Dew T and composition

% P = 1atm = 101.325 kPa
P = 101.325;
epsilon = 1e-6;

% Initialize 
phi (1:n) = 1;
gamma(1:n) = 1;
y_new = [0.17, 0.23, 0.29, 0.31];
n = length(y_new);
x_new = zeros(1,n);
 
% Calculate Tsat
T_sat = Cal_T(P,n);

%Calc T
T1 = y_new*T_sat';

% Calc Psat 
P_sat = Cal_P(T1,n);

%Toluene is used as species j

Pj_sat =P*((y_new * P_sat(1) .* phi ./gamma)*(P_sat.^(-1))');

%Calc T
T2= (B/(A - log(Pj_sat)))-C;

%define delta for termination criteria
delta_T= abs(T2 - T1);

while delta_T >= eps
    T1 = T2; %set T_1 to old T for delta calc later
    P_sat = Cal_P(T1,n); %recalculate P_sat for new T
    
    %calculate liquid mixture composition
    for i=1:n
        x_new(i)=y_new(i)*phi(i)*P/(gamma(i)*P_sat(i));
    end
    
    normalise = sum(x_new);
    %normalize x if neccessary
    if normalise ~= 1
        x_new = x_new*(1/normalise);
    end
    
    %calc new Pj_sat
    Pj_sat = P*((y_new * P_sat(1) .* phi ./gamma)*(P_sat.^(-1))');
    
    %calc new T
    T2= (B/(A-log(Pj_sat)))-C;
    
    %calc new delta
    delta_T= abs(T2 - T1);
    
end

SOL_q2b_T = T2;
SOL_q2b_x = x_new;



%% Problem 3
% - Assume non-ideal behavior, fully liquid, calculate Bubble T and
%  composition in the gas phase.
% - Assume non-ideal behavior, fully gas and calculate Dew T and 
% composition in the liquid phase


%all values from table B.1 in 7th ed
           %tol    Benze  Cyclo  Chloro
omega_pure  = [ 0.262, 0.210, 0.210, 0.250];
Tc_pure = [ 591.8, 562.2, 553.6, 632.4]; %in kelvin
Zc_pure = [ 0.264, 0.271, 0.273, 0.265];
Vc_pure = [ 316.0, 259.0, 308.0, 308.0]; %in cm^3/mol
x_new = [0.17, 0.23, 0.29, 0.31];
n = length(x_new);

% ===========Bubble point temperature=============
%initialize phi and y
phi_bubl(1:n)= 1; 
y_bubl = zeros(1,n);


%calculate Tsat
Tsat_bubl= Cal_T(P,n);

% calc T
T1_bubl = x_new * Tsat_bubl';

% calculate Psat 
P_sat = Cal_P(T1_bubl,n);

% calculate gamma
gamma_bubl3 = UNIFAC(T1_bubl + 273.15,x_new);

% Species j: Toluene
Pj_sat = P/sum((x_new.*gamma_bubl3.*P_sat)./(P_sat(1)*phi_bubl));
T2_bubl = (B/(A-log(Pj_sat)))-C;

deltaT = abs(T2_bubl-T1_bubl);

while (deltaT >= epsilon)
    T1_bubl = T2_bubl;
    P_sat = Cal_P(T2_bubl,n);
    y_bubl=(x_new.*gamma_bubl3.*P_sat)./(P * phi_bubl);
    
    [phi_bubl] = bigPhi_calc(P,P_sat,T2_bubl+273.15,y_bubl,omega_pure,Tc_pure,Zc_pure,Vc_pure);
    gamma_bubl3=UNIFAC(T2_bubl+273.15,x_new);
    
    Pj_sat=P/sum((x_new.*gamma_bubl3.*P_sat)./(P_sat(1)*phi_bubl));
    
    T2_bubl=(B/(A-log(Pj_sat)))-C;
    
    deltaT=abs(T2_bubl - T1_bubl);
  
end

SOL_q3_bublT = T2_bubl;
SOL_q3_y = y_bubl;


% =========== Dew point temperature=============
y_dew = [0.17, 0.23, 0.29, 0.31];
eta = 1e-6;
n=length(y_dew);
epsilon = 1e-6;
%initialize variable matrices
bigPhi_dew(1:n)= 1;
gamma_dew(1:n) = 1; 
%x_dew = zeros(1,n);

Tdew_sat = Cal_T(P,n);

T1_dew = y_dew*Tdew_sat';

%calc P_sat for each species
Pdew_sat = Cal_P(T1_dew,n);

% Species j: Toluene
Pjsat_dew = P*sum((y_dew.*bigPhi_dew.*Pdew_sat(1))./(gamma_dew.*Pdew_sat));


T2_dew =(B/(A-log(Pjsat_dew)))-C; % in C


Pdew_sat = Cal_P(T2_dew,n);
[bigPhi_dew]= bigPhi_calc(P, Pdew_sat,T2_dew + 273.15,y_dew,omega_pure,Tc_pure,Zc_pure,Vc_pure);
x_dew=(y_dew.*bigPhi_dew*P)./(gamma_dew.*Pdew_sat);
new_gamma = UNIFAC(T2_dew + 273.15,x_dew);

Pjsat_dew=P*sum((y_dew.*bigPhi_dew.*Pdew_sat(1))./(gamma_dew.*Pdew_sat));

T2_dew =(B/(A-log(Pjsat_dew)))- C;
deltaT_dew = abs(T2_dew - T1_dew);


while deltaT_dew > epsilon 
    T1_dew = T2_dew;
    Pdew_sat = Cal_P(T2_dew,n);
    [bigPhi_dew]= bigPhi_calc(P, Pdew_sat,T2_dew + 273.15,y_dew,omega_pure,Tc_pure,Zc_pure,Vc_pure);
    
  
    delta_gamma = 1;
    while delta_gamma > eta
        
        % Calculate xdew
        x_dew=(y_dew.*bigPhi_dew*P)./(gamma_dew.*Pdew_sat);
        
        %normalize x if neccessary
        if sum(x_dew)~= 1
            x_dew = x_dew/sum(x_dew);
        end
        new_gamma = UNIFAC(T2_dew + 273.15,x_dew);
        delta_gamma = abs(new_gamma-gamma_dew);
        
        gamma_dew = new_gamma;
    end
   % fprintf("gamma: \n");
    %disp(new_gamma);
    Pjsat_dew=P*sum((y_dew.*bigPhi_dew.*Pdew_sat(1))./(gamma_dew.*Pdew_sat));

    T2_dew =(B/(A-log(Pjsat_dew)))- C;
    deltaT_dew = abs(T2_dew - T1_dew);
    
end

SOL_q3_dewT = T2_dew;
SOL_q3_x = x_dew;

%% Problem 4: Vapor and Liquid composition
% Assume system at 373K and 1 atm, calculate the mole fraction of the
% liquid and vapor phase, and their compositions


%all values from table B.1 in 7th ed textbook
           %tol    Benze  Cyclo  Chloro
omega_pure  = [0.262, 0.210, 0.210, 0.250];
Tc_pure = [ 591.8, 562.2, 553.6, 632.4]; %in kelvin
Zc_pure = [ 0.264, 0.271, 0.273, 0.265];
Vc_pure = [ 316.0, 259.0, 308.0, 308.0]; %in cm^3/mol

x_new = [0.17, 0.23, 0.29, 0.31];

T = 373; % Kelvin
z = [0.17, 0.23, 0.29, 0.31];
n = length(z);
P_sat4 = Cal_P(T - 273, n);
epsilon = 1e-6;
eta = 1e-6;
y_dew4 = z;
x_dew4 = zeros(1,n);

bigPhi_dew4(1:n)= 1; 
gamma_dew4(1:n)= 1; 

%calc Dew P
P_dew4 = 1/(y_dew4*(1./P_sat4)');

deltaP_dew = abs(P_dew4 - P);

while deltaP_dew > epsilon
    P_dew_new = P_dew4; 
    [bigPhi_dew4]= bigPhi_calc(P_dew4, P_sat4, T, y_dew4,omega_pure,Tc_pure,Zc_pure,Vc_pure);
    delta_gamma = 1;
    while delta_gamma > eta 
        
        x_dew4 = y_dew4 .* bigPhi_dew4* P_dew4 ./(gamma_dew4.*P_sat4); 
       
        if sum(x_dew4)~= 1
            x_dew4 = x_dew4 / sum(x_dew4);
        end
        
        % Calculate gamma
        gamma_dew_new4 = UNIFAC(T,x_dew4);
        
        delta_gamma = abs(gamma_dew4 - gamma_dew_new4);
        gamma_dew4 = gamma_dew_new4;
    end
    
    %Calculate P_dew
    
    term = y_dew4.*bigPhi_dew4./(gamma_dew4.*P_sat4);     
    P_dew4 = 1/sum(term);
    deltaP_dew = abs(P_dew4 - P_dew_new);
end
   
%calc P bubble
x_bubl4 = z;
bigPhi_bubl4(1:n)= 1;

gamma_bubl4 = UNIFAC(T,x_bubl4);
P_bubl4 = sum(x_bubl4.*gamma_bubl4.*P_sat4/bigPhi_bubl4);

delta_Pbub = abs(P_bubl4 - P);

while delta_Pbub > epsilon
    %Calculate y bubl
     y_bubl4 = x_bubl4 .*gamma_bubl4.*P_sat4./(bigPhi_bubl4.*P_bubl4);
     
     %Calculate big Phi
    [bigPhi_bubl4] = bigPhi_calc(P_bubl4,P_sat4,T,y_bubl4,omega_pure,Tc_pure,Zc_pure,Vc_pure);

    term = x_bubl4.*gamma_bubl4.*P_sat4./bigPhi_bubl4;
    
    P_bubl_new = sum(term);
    delta_Pbub = abs(P_bubl_new - P_bubl4);
    P_bubl4 = P_bubl_new;
    
end

if (P < P_dew4 || P > P_bubl4)
    error('Pressure is not in the range.Stop.\n');
end

%estimate V, Phi, and gamma

V_4 = 1 + (P_dew4 - P)/(P_bubl4 - P_dew4);
gamma_4 = gamma_dew4 + (V_4 *(gamma_bubl4 - gamma_dew4));
bigPhi_4 = bigPhi_dew4 + (V_4 *(bigPhi_bubl4 - bigPhi_dew4));

delta_x = norm(x_dew4 - x_bubl4,2); 
delta_y = norm(y_bubl4 - y_dew4,2);
delta_V = 1;

x_4 = x_bubl4;
y_4 = y_dew4;

while (delta_x > epsilon || delta_y > epsilon || delta_V > epsilon)
    K_4 = gamma_4.* P_sat4 ./ bigPhi_4 / P;
    approx = 1;
    while (approx > epsilon)  
       
        F_4 = sum((z.*(K_4 -1 ))./(1 + V_4*(K_4 - 1)));
        dF_4 = -sum((z.*((K_4 - 1).^2))./((1 + V_4*(K_4 - 1)).^2));
        V_new4 = V_4 - (F_4/dF_4);
        approx = abs(F_4/dF_4);
        V_4 = V_new4;
    end
    %Calculate newx, new y
    x_new = z./(1 + V_4*(K_4-1));
    y_new = K_4 .* x_new;
    
    %Calculate gamma
    gamma_4 = UNIFAC(T,x_new);
    
    %Calculate bigPhi
    [bigPhi_4]= bigPhi_calc(P, P_sat4, T, y_new,omega_pure,Tc_pure,Zc_pure,Vc_pure);
    
    % Calculate delta x, delta y, and delta V
    delta_x = norm(x_new - x_4,2);
    x_4 = x_new;
    delta_y = norm(y_new - y_4,2);
    y_4 = y_new;
    
    delta_V = abs(V_new4 - V_4);
    
end

SOL_q4_y = y_new;
SOL_q4_x = x_new;
SOL_q4_V = V_4;
SOL_q4_L = 1 - V_4;

%% Problem 5:
% Calculate T and P given such that the liquid phase composition of toluene y1=0.10

eps = 1e-6; %termination criteria

%Antoine's Eqn Coeff for Toluene


x=[0.17, 0.23, 0.29, 0.31];
y=[0.10, 0, 0, 0];
n=length(x);

T = 400; %K<--guess
bigPhi(1:n)=1; 

P_sat = Cal_P(T - 273.15,n);
gamma = UNIFAC(T,x);
P_5 = sum((x.*gamma.*P_sat)./bigPhi);


delta_P = 1;
   
while delta_P > eps
   T_new =(B/(A-log(P_sat(1))))-C; %(in C degree)
   delta_T = abs(T_new -T); 
    
    while delta_T > eps
        % Calculate y
       y(2:n)=(x(2:n).*gamma(2:n).*P_sat(2:n))./(bigPhi(2:n)*P_5);
       
       % Calculate Psat of toluene
       P_sat(1)=y(1)*bigPhi(1)*P_5/(x(1)*gamma(1));
       
       % Calc T_new
       T_new =(B/(A-log(P_sat(1))))-C; %(in C degree)
       % Calc Pisat
       P_sat = Cal_P(T_new,n);
       
       %Calc gamma and bigPhi
       gamma = UNIFAC(T_new+ 273.15, x);     
       [bigPhi] = bigPhi_calc(P_5,P_sat,T_new+273.15,y,omega_pure,Tc_pure,Zc_pure,Vc_pure);
       
       % Calc delta T and then update T
       delta_T = abs(T_new+273.15 - T);
       T = T_new + 273.15;
    
    end
    
    % Use gamma, Psat, bigPhi at that T to calcualte P
    P_new = sum((x.*gamma.*P_sat)./bigPhi);
    delta_P = abs(P_5-P_new);
    P_5 = P_new;
    
end

SOL_q5_T = T - 273.15; % in C
SOL_q5_P= P_5 * 1e3; % in Pa
SOL_q5_y = y;


%% Functions

% Caluclate Psat and Tsat using Antoine equation. 

function P_sat = Cal_P(T,n)
% Input: Temperature(T in C degree) and number of species (n)
% Output: P_sat(kPa) for each species
% Method: Antoine equation

% Antoine's Coefficients for 4 species:
% toluene, benzene, cyclohexane, chlorobenzen

A = [13.9320, 13.7819, 13.6568, 13.8635];
B = [3056.96, 2726.81, 2723.44, 3174.78];
C = [217.625, 217.572, 220.618, 211.700];
% Initialize P_sat array
P_sat = zeros(1,n);

% Calculate P_sat
    for i=1:n
       P_sat(i)= exp(A(i)-(B(i)/(T+C(i))));
    end

end


function T_sat = Cal_T(P,n)
% Input: Pressure(P) and number of species (n)
% Output: T_sat(K) for each species
% Method: Antoine equation

% Antoine's Coefficients for 4 species:
% toluene, benzene, cyclohexane, chlorobenzen

A = [13.9320, 13.7819, 13.6568, 13.8635];
B = [3056.96, 2726.81, 2723.44, 3174.78];
C = [217.625, 217.572, 220.618, 211.700];

% Initialize T_sat array
T_sat = zeros(1,n);

% Calculate T_sat
    for i=1:n
        T_sat(i)= B(i)/(A(i)-log(P)) - C(i);
    end

end


% Calculate gamma using UNIFAC

function gamma = UNIFAC(T, x)
% Input: Temp(T in K), mole fraction (x)
% Output: A array of gamma for each species based on temp and
% liquid compositions
% Method: UNIFAC

% UNIFAC group parameters
%Species: (vk)*subgroup
% 1- Toluene:       5*ACH(3), 1*ACCH3(4)
% 2- Benzene:       6*ACH(3)
% 3- Cyclohexane:   6*CH2(1)
% 4- Chlorobenzene: 5*ACH(3), 1*ACCL(25)

n = length (x);
%Rk and Qk for subgroups
Rk_CH2  = 0.6744; 
Qk_CH2  = 0.540;

Rk_ACH  = 0.5313;
Qk_ACH  = 0.400;

Rk_ACCH3  = 1.2663; 
Qk_ACCH3  = 0.968;

Rk_ACCl = 1.1562; 
Qk_ACCl = 0.844;

% Calculate r and q for species

r =[5*Rk_ACH + Rk_ACCH3, 6*Rk_ACH, 6*Rk_CH2, 5*Rk_ACH + Rk_ACCl];

q =[5*Qk_ACH + Qk_ACCH3, 6*Qk_ACH, 6*Qk_CH2, 5*Qk_ACH + Qk_ACCl];



% calculate eki=vk(i)Qk/qi
%species     1              2              3                4
e_k = [              0,             0, 6*Qk_CH2/q(3),               0;  
         5*Qk_ACH/q(1), 6*Qk_ACH/q(2),             0,   5*Qk_ACH/q(4);  
       1*Qk_ACCH3/q(1),             0,             0,               0;  
                     0,             0,             0, 1*Qk_ACCl/q(4)];  
                 

             
% J = ri/sum(rj*xj) (Eq G.10)
J = r/(x*r'); 

% L = qi/sum(qj*xj) (Eq G.11)
L = q/(x*q');


%initialize amk matrix with values from the UNIFAC website

%subgroup  1      3        4         25
a_mk =  [     0,   61.13,  76.50,   11.44;  %  1
        -11.12,       0, 167.00,  187.00;  %  3
        -69.70, -146.80,      0, -211.00;  %  4
        106.80,  -97.27, 402.50,       0]; % 25

a_mk = a_mk';
% calculate tau=exp (-amk/T) (Eq G. 21)
tau = exp(-a_mk/T);

% calculate beta_ik = summ(e_m*tau_mk) (Eq G. 18)
beta_ik = tau*e_k;
beta_ik = beta_ik';


% calculate theta (Eq G.19)
theta_num = zeros(n);
for k=1:n
   
   theta_num(:,k)=x.*q.*e_k(k,:);
    
end
theta = sum(theta_num)/ sum(x.*q);


%calc sk
s_k = theta(:,1:n)*tau';

%calculate gamma_C
gamma_C = exp(1-J+log(J)-(5*q.*(1-(J./L)+log(J./L)))); %G.13)


%calculate gamma_R
gamma_R = zeros(1,n);  

e_k = e_k';
 for i=1:n 
    term = 1 - sum(theta .*beta_ik(i,:)./s_k - e_k(i,:).*log(beta_ik(i,:)./s_k), 'all');
    gamma_R(i)= exp(q(i)*term);
    
 
 end



%calculate gamma
gamma = gamma_C.*gamma_R;


end

% Calculate bigPhi

function [big_phi] = bigPhi_calc(P,Psat,T,y,omega_pure,Tc_pure,Zc_pure,Vc_pure)
% Calculate big Phi
% T in Kelvin


n = length(y);
R = 8314; % ideal gas constant with V in cm^3


%initialize matrices for mixture variables
T_c = zeros(n);
omega = zeros(n); 
Z_c = zeros(n); 
delta = zeros(n); 
V_c = zeros(n);
B_1 = zeros(n); 
B_0 = zeros(n); 
B = zeros(n); 



for i=1:n
    for j=1:n
        T_c(i,j) = (Tc_pure(i)*Tc_pure(j))^0.5;
        omega(i,j)  = (omega_pure(i)+omega_pure(j))/2;
        Z_c(i,j) = (Zc_pure(i)+Zc_pure(j))/2;
        V_c(i,j) = ((Vc_pure(i)^(1/3)+Vc_pure(j)^(1/3))/2)^3;
        B_0(i,j) = 0.083 - (0.422/((T/T_c(i,j))^(1.6)));
        B_1(i,j) = 0.139 - (0.172/((T/T_c(i,j))^(4.2)));
        B(i,j)  = V_c(i,j)*(B_0(i,j)+(omega(i,j)*B_1(i,j)))/Z_c(i,j);
        
        
    end
end

% calculate delta
for i=1:n
    for j=1:n
        delta(i,j)= 2*B(i,j)-B(i,i)-B(j,j);
    end
end

%calc Phi

term1 = zeros(1,n);
term2 = zeros(1,n);
big_phi = zeros(1,n);



% Calculate bigPhi
for i=1:n 
   for j=1:n
        for k=1:n
            term1(k)=(y(j)*y(k)*((2*delta(i,j))-delta(j,k)));
        end
   term2(j)=sum(term1);
   end
   big_phi(i) = exp((B(i,i)*(P - Psat(i))+ 0.5*P*sum(term2))/(R*T));
end
    
   

end
