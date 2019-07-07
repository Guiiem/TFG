clc
close all
clear variables

addpath ..\01_Non_iterative_development\;
tic

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Generation of the initial population
[Chord,Twist] = Creator(Data);

%%Compute the AEP for each individual
for i=1:Data.no_indiv
    
    %%BEMT calculation for a range of lambda
    lambda = 1:0.25:10;
    for j=1:length(lambda)
        Data.lambda_power = lambda(j);
        [Cp(j), Cq(j)] = power_calc_opt(Data,Chord(i,:),Twist(i,:),0);
    end

    %%Torque and Power vs wind speed as a function of Lambda
    Lambda_Results = Curves(Cp,Cq,Data,lambda,Plots);

    %%Find the intersection point between results and generator curve
    Operation = FindOp(Lambda_Results,Data,lambda,Plots);

    %%With the obtained results, compute the Annual Energy Production
    AEP(i) = AnnualEnergyProduction(Operation,Data,Plots);
end

simulation_time = toc

