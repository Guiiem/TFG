clc
clear variables
close all

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Creation of the chord and twist distribution
[Chord,Twist] = ChordTwist(Data);

%%BEMT calculation for a range of lambda
lambda = 1:0.25:10;
for i=1:length(lambda)
    Data.lambda_power = lambda(i);
    [Cp(i), Cq(i)] = power_calc_opt(Data,Chord,Twist,0);
end

%%Torque and Power vs wind speed as a function of Lambda
Results = Curves(Cp,Cq,Data,lambda,Plots);

%%Find the intersection point between results and generator curve
Operation = FindOp(Results,Data,lambda,Plots);

%%With the obtained results, compute the Annual Energy Production
AEP = AnnualEnergyProduction(Operation,Data,Plots);







    

