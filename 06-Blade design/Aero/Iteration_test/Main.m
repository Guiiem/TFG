% Wind Turbines Desgin Project 
% Guillem Vergés i Plaza
% Main

clearvars;
close all;
%% 1. Data input

Nb = 3; %Number of blades
A = 1; %Last DNI numerical digit
D = 1.2*(A/2 + 1.5); %Rotor diameter [m]
c = 0.08*(A/2 + 1); %Blade chord [m]
Umax = 25; %Maximum wind speed in operation [m/s]
Rho = 1.225; %Air density [kg/m^3]
Ad = pi*(D/2)^2; %Swept area [m^2]

% Discretization

Mu = [0.4 0.6 0.8 0.95]; %Position of each node [m/m]
r = D/2*Mu; %Position of each node [m]
Delta_r = D/2*[0.2 0.2 0.2 0.1]; %Length of each node [m]
tau = [6-6.67*0.4 6-6.67*0.6 6-6.67*0.8 2]; %Twist for each node [º]
Sigma_r = (Nb*c)./(2*pi*r); %Chord solidity 

Cp_opt = 0;
Delta = 1e-5; %Convergence criteria




%% 2. Maximum power coefficient 
%We need to compute the power coefficient for each 
%combination of tip speed ratio and pitch angle

%Tip speed ratios and pitch angles that will be studied
theta = 0:2:8;
lambda = 5:2:15;

%Creation of the matrices 
CP = zeros (length(lambda), length(theta));
CT = zeros (length(lambda), length(theta));

%Iteration for each tip speed ratio    
for j=1:length(lambda);

    %Iteration for each pitch angle
    for k=1:length(theta);
        %We use a function that iterates with those given values to
        %find the coefficient of power and thrust
        [CP(j,k),CT(j,k)] = PowerCoef(theta(k), tau, lambda(j), Mu, Sigma_r, Delta, r, Delta_r, D);

        %Evaluate if this is a new maximum
        if (CP(j,k) > Cp_opt)
            Cp_opt = CP(j,k);
            lambda_opt = lambda(j);
            theta_opt = theta(k);
        end        
    end
end

Cp_opt = 0.5316;
lambda_opt = 7.83;
theta_opt = 0.33;

%Plot the coefficient of power vs theta for different lambdas

figure('Name','Cp')
for i=1:length(lambda)
    plot(theta,CP(i,:));
    ylim([0 0.6]);
    legendInfo{i}=['\lambda: ' num2str(lambda(i))];
    hold on
end
xlabel('Angle \theta [º]');
ylabel('Coefficient of power C_p');
grid on;
grid minor;
title('Coefficient of power as function of \theta and \lambda');
legend(legendInfo);

%% 3. Noise limitation and rotational speed

w_max = 24/(A/2 + 1.5); %Maximum rotational speed [Hz=1/s]
omega_max = w_max*2*pi; %Maximum rotational speed [Rad/s]
U_max = omega_max*(D/2)/lambda_opt; %Wind speed for omega max

U_inf = 0.01:0.2:25; %Range of wind speed 
omega = length(U_inf); %Range of rotor speed
lambda_real = length(U_inf); %Range of actual lambda 


for i=1:length(U_inf)
    omega(i) = lambda_opt*U_inf(i)/(D/2); %Rotational speed [Rad/s]
    if omega(i) > omega_max
        omega(i) = omega_max;
    end
    lambda_real(i) = omega(i)*(D/2)/U_inf(i);
end

figure ('Name','Rotational speed')
yyaxis left
plot(U_inf,omega*60/(2*pi)) %For making the plot in RPM
ylabel('Rotational speed \omega [rpm]')
hold on
yyaxis right
plot(U_inf,lambda_real)
ylabel('Tip speed ratio')
title('Rotational speed with noise limitations')
h = line([U_max U_max], [0 9],'Color','k','LineStyle','--');
grid on;
grid minor;
xlabel('Wind speed [m/s]')

%% 4. Operation curves

NP = 2.5*A*1000; %Nominal power [W]
P = length(U_inf); %Power at each wind speed [W]
Cp_real = length(U_inf); %Real Cp of the wind turbine
theta_real = length(U_inf); %Real angle theta
Ct = length(U_inf); %Coefficient of trhust vector
T = length(U_inf); %Thrust vector

%PROVA
for i=1:length(U_inf)
        [CP,CT] = PowerCoef(theta_opt, tau, lambda_real(i), Mu, Sigma_r, Delta, r, Delta_r, D);
        Cpprova(i) = CP;
end
figure()
plot(U_inf,Cpprova)
xlabel('Wind speed amb la seva lambda corresponent')
ylabel('Cp')


 for i=1:length(U_inf)
    %Evaluate the axial induction factor and the coefficient of power 
    %taking into account the noise limitations
    [CP,CT] = PowerCoef(theta_opt, tau, lambda_real(i), Mu, Sigma_r, Delta, r, Delta_r, D);
    P(i) = CP*0.5*Rho*Ad*U_inf(i)^3; %Calculation of power
    T(i) = CT*0.5*Rho*Ad*U_inf(i)^2; %Calculation of thrust
    theta_real(i) = theta_opt; 
    Cp_real(i) = CP;
    Ct(i) = CT;
     
    if P(i)> NP || P(i)~=P(i)
        if( P(i)~=P(i))
        power = P(i)
        winspit = U_inf(i)
        end
        P(i) = NP; %If the power reaches nominal power, we will vary the pitch
        %to keep NP constant
        
        %We will have a lower Cp, as theta should have increased
        Cp_real(i) = P(i)/(0.5*Rho*Ad*U_inf(i)^3); 
        
        %Finding theta for this Cp (ALGUNA MANERA MILLOR DE FER-HO?)
        Theta_find= 0:1:20; %Range of theta 
        dif = 10; 
        
        for j=1:length(Theta_find)
            [CP,CT] = PowerCoef(Theta_find(j), tau, lambda_real(i), Mu, Sigma_r, Delta, r, Delta_r, D);
            if abs(CP-Cp_real(i)) < dif
                theta_real(i) = Theta_find (j);
                dif = abs(CP-Cp_real(i));   
                
                %Actualizing the trhust
                T(i) = CT*0.5*Rho*Ad*U_inf(i)^2; 
                Ct(i) = CT;             
            end           
        end        
    end  
 end

 %Power
 figure('Name','Power and thrust')
 plot(U_inf,P,'Color',[0, 0.4470, 0.7410])
 hold on
 plot(U_inf,T,'LineStyle','--','Color',[0, 0.4470, 0.7410])
 ax = gca;
 ax.YColorMode = 'manual';
 ax.YColor = [0, 0.4470, 0.7410]	;
 addaxis(U_inf,Cp_real,[0 1],'Color',[0.8500, 0.3250, 0.0980]);
 hold on
 addaxisplot(U_inf,Ct,2,'LineStyle','--','Color',[0.8500, 0.3250, 0.0980])
 hold on
 addaxis(U_inf,theta_real,[0 30],'Color',[0.4660, 0.6740, 0.188],'LineStyle','-')
 grid on
 grid minor
 hold on
 ylabel('Power and Thrust [N]');
 addaxislabel(2,'Power and Thrust Coefficients');
 addaxislabel(3,'Theta angle \theta [º]');
 legend('P','T','C_P','C_T','\theta')
 xlabel('Wind speed [m/s]')
 title('Power curve characteristics')
 
 
%% 5. Annual energy production (AEP)

%Weibull distribution parameters
A = 7; %[m/s]
k = 1.9;

Prob = k/A .* (U_inf/A).^(k-1) .* exp(-(U_inf/A).^k); %Compute the probabilty of the Weibull distribution
Hours = Prob*24*365; 

figure('Name','AEP')
yyaxis left
plot(U_inf,Prob)
ylabel('Probability')
xlabel('Wind speed [m/s]')
yyaxis right
plot(U_inf,P)
ylabel('Power output [kW]')
grid on;
grid minor;
title('Power and probability for each wind speed')

Energy_produced = Prob.*P*24*365; %Energy produced for each wind speed in one year [W*h]



%Now we can integrate numerically the energy produced
deltaU=U_inf(2)-U_inf(1); %Diferential of wind speed
deltaE = deltaU*Energy_produced; %Contribution of each diferential of energy
sumE = length(U_inf)
for i=1:length(U_inf)
    if i==1
        sumE(i) = deltaE(i);
    end
    if i~=1
    sumE(i) = sumE(i-1) + deltaE(i);
        end
end
AEP = sum(deltaE)

AEP = trapz(U_inf,Energy_produced)

figure('Name','AEP2')
yyaxis left
plot(U_inf,Energy_produced/1000000)
ylabel('Energy produced [MW h]')
yyaxis right
plot(U_inf,sumE/1000000) 
ylabel('Acumulated energy production [MW h]')
grid on
grid minor
xlabel('Wind speed [m/s]')
title('Annual Energy Production')













