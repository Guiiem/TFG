function [AEP, Data, Operation, Power_Curve] = BigSolver (Data, Chord, Twist, Plots)

%%BEMT calculation for a range of lambda
lambda = Data.lambda_range; 
Power_Curve = struct; %Create struct
for j=1:length(lambda)
    Data.lambda_power = lambda(j);
    %[Cp(j), Cq(j)] = power_calc_opt(Data, Chord, Twist,0); %Book code
    Power_Curve = PowerCoef_v2_WithTipLosses(Power_Curve,Data,Chord,Twist,j); %Own code
end

%%Torque and Power vs wind speed as a function of Lambda
Lambda_Results = Curves(Power_Curve.CP, Power_Curve.CQ, Power_Curve.CT, Data,lambda,Plots);

%%Find the intersection point between results and generator curve
Operation = FindOp(Lambda_Results,Data,lambda,Plots,Power_Curve);

%%With the obtained results, compute the Annual Energy Production
AEP = AnnualEnergyProduction(Operation,Data,Plots);

if Plots.CT == true
    figure('Name','Chord and twist distribution')
    subplot(121)
    plot(Data.rad*Data.r_tip,Chord)
    xlabel('Radius [m]')
    ylabel('Chord [m]')
    
    subplot(122)
    plot(Data.rad*Data.r_tip,Twist)
    xlabel('Radius [m]')
    ylabel('Twist [º]')
end

