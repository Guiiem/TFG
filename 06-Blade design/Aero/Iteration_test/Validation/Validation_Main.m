clc
close all
clear variables

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Extrapolate the 360º polar curves 
Data = Viterna_Extrapolation(Data);

%%Select which plots you want to be printed or not
Plots = PlotSelectionFinal;

%%Generation of the twist and chord distribution to be tested
Blade = Creator_Validation(Data);

%Power_Curve = struct; %Create struct
Power_Curve = struct;

%Evaluate Cp_Lambda
for i=1:length(Data.lambda_range)
    Data.lambda_power = Data.lambda_range(i);
    Power_Curve = PowerCoef_v2_WithTipLosses(Power_Curve,Data,Blade.Chord(1,:),Blade.Twist(1,:),i); %Own code with tip losses
end

lambda = Data.lambda_range;
Lambda_Results = Curves(Power_Curve.CP, Power_Curve.CQ, Power_Curve.CT, Data,lambda,Plots);

