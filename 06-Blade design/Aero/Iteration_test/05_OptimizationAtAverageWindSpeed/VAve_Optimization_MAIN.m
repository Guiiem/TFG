clc
% close all
clear variables

addpath ..\01_Non_iterative_development\;
addpath ..\02_Genetic_algorithm_development\;
addpath ..\03_PowerCalcTest\;
addpath ..\04_TradeOff_StartingAEP\;

%Optimization wind speed
U_opt = [3 5 7 9 10 11]; %[m/s]


%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Extrapolate the 360� polar curves 
Data = Viterna_Extrapolation(Data);

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Generation of the twist and chord distribution to be tested
Blade = Creator_VAve_Optimization(Data);

%Power_Curve = struct; %Create struct
Power_Curve_TL = struct;


%Compute the Cp and power for each blade geometry
for i=1:length(Blade.Lambda)
    Data.lambda_power = Blade.Lambda(i);
    %Power_Curve = PowerCoef_v2_WithOUTTipLosses(Power_Curve,Data,Blade.Chord(i,:),Blade.Twist(i,:),i); %Own code
    Power_Curve_TL = PowerCoef_v2_WithTipLosses(Power_Curve_TL,Data,Blade.Chord(i,:),Blade.Twist(i,:),i); %Own code
    %[cp(i),torque(i),aoa(i,:),a(i,:)] = power_calc_opt_VAveOptimization(Data,Blade.Chord(i,:),Blade.Twist(i,:), 0);
    %P(i) = Power_Curve.CP(i)*0.5*1.225*Data.r_tip^2*pi*U_opt^3;
    
    for j = 1:length(U_opt)
        P_TL(i,j) =Power_Curve_TL.CP(i)*0.5*1.225*Data.r_tip^2*pi*U_opt(j)^3;
        Omega(i,j) = U_opt(j)*Blade.Lambda(i)/Data.r_tip *60/2/pi; %Rotor speed in rpm
        %P(i) = cp(i)*0.5*1.225*Data.r_tip^2*pi*U_opt^3;
    end
end

figure()
plot(Omega,P_TL)
hold on
plot(Data.gen_omega, Data.gen_power)
legend('Uopt = 3 m/s','Uopt = 5 m/s','Uopt = 7 m/s','Uopt = 9 m/s','Uopt = 10 m/s','Uopt = 11 m/s','Generator curve')
grid on
xlabel('Rotor speed [rpm]')
ylabel('Power [W]')

%Find intersection lambda








