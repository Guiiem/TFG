clc
close all
clear variables

%Read QBlade Data
Read_QBlade_Data;
Read_ClArray;
Read_CdArray;

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;


%%Extrapolate the 360º polar curves 
%Data = Viterna_Extrapolation(Data);

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Generation of the twist and chord distribution to be tested
Blade = Creator_Validation(Data);


Data.cl_array = polarcl(:,2);
Data.cd_array = polarcd(:,2);
Data.alpha = polarcl(:,1);



%Power_Curve = struct; %Create struct

Hansen_QB = struct;
Hansen_QB_TipLoss = struct;
Hansen_QB_TipRootLoss = struct;


%Evaluate Cp_Lambda
for i=1:length(Data.lambda_range)
    Data.lambda_power = Data.lambda_range(i);
    Options.TipLoss = 0;
    Options.RootLoss = 0;
    Hansen_QB = PowerCoef_Hansen_QBlade_Total(Options,Hansen_QB,Data,Blade.Chord_Boundary(1,:),Blade.Twist_Boundary(1,:),i);
    Options.TipLoss = 1;
    Hansen_QB_TipLoss = PowerCoef_Hansen_QBlade_Total(Options,Hansen_QB_TipLoss,Data,Blade.Chord_Boundary(1,:),Blade.Twist_Boundary(1,:),i);
    Options.RootLoss = 1;
    Hansen_QB_TipRootLoss = PowerCoef_Hansen_QBlade_Total(Options,Hansen_QB_TipRootLoss,Data,Blade.Chord_Boundary(1,:),Blade.Twist_Boundary(1,:),i);

    
end

lambda = Data.lambda_range;
%Lambda_Results_OwnCode_NoTipLosses = Curves(Power_Curve_OwnCurve_NoTipLosses.CP, Power_Curve_OwnCurve_NoTipLosses.CQ, Power_Curve_OwnCurve_NoTipLosses.CT, Data,lambda,Plots);
%Lambda_Results_OwnCode_TipLosses = Curves(Power_Curve_OwnCurve_TipLosses.CP, Power_Curve_OwnCurve_TipLosses.CQ, Power_Curve_OwnCurve_TipLosses.CT, Data,lambda,Plots);
%Lambda_Results_BookCode = Curves(Book.Cp, Book.Cq, Book.Ct, Data, lambda, Plots);

figure()
plot(Hansen_QB.lambda,Hansen_QB.Cp)
hold on
plot(Hansen_QB_TipLoss.lambda,Hansen_QB_TipLoss.Cp)
plot(Hansen_QB_TipLoss.lambda,Hansen_QB_TipRootLoss.Cp)
plot(Cplambdavalidation(:,1),Cplambdavalidation(:,2))
plot(Cplambdavalidation(:,1),Cplambdavalidation(:,4))
plot(Cplambdavalidation(:,1),Cplambdavalidation(:,6))
legend('Own code','Own code - Tip Loss','Own code - Tip + Root Loss','Qblade','Qblade - Tip Loss' ,'Qblade - Tip + Root Loss')
xlim([0.5 3.5])

figure()
hold on
%plot(Hansen_TipLoss_QB.rad
%plot(Hansen.lambda,Power_Curve_OwnCurve_NoTipLosses.CP)
%plot(Hansen.lambda,Power_Curve_OwnCurve_TipLosses.CP)
plot(Hansen_QB.rad,Hansen_QB.a(10,:));
plot(Hansen_QB.rad,Hansen_QB_TipLoss.a(10,:));
plot(Hansen_QB.rad,Hansen_QB_TipRootLoss.a(10,:));
legend('Hansen code','Hansen code - Tip Loss','Hansen/Qblade code - Tip Loss')

