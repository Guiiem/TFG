clc
close all
clear variables

addpath ..\01_Non_iterative_development\;
addpath ..\02_Genetic_algorithm_development\;
addpath ..\03_PowerCalcTest\;
addpath ..\05_OptimizationAtAverageWindSpeed\;

tic

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Extrapolate the 360º polar curves 
Data = Viterna_Extrapolation(Data);

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Generation of the twist and chord distribution to be tested
Blade = Creator_TradeOff_v3(Data);

%Compute the AEP for each individual
for i=1:Data.no_indiv
   [Blade.AEP(i), Data, Operation(i),~] = BigSolver (Data, Blade.Chord(i,:), Blade.Twist(i,:), Plots);
   Blade.J(i) = inertia(Blade.Chord(i,:), Blade.Twist(i,:), Data);
   Blade.Ts(i) = Start_Calc_Opt_v2(Data, Blade.Chord(i,:), Blade.Twist(i,:), Blade.J(i), Operation(i));
end

%Select the best option
Ind = TradeOff_Selection(Blade, Data, Operation);

%Plot the results
Plots = PlotSelectionFinal; %Now we want to see the final plots

[AEP, ~, Operation_selected,Power_Curve] = BigSolver(Data, Blade.Chord(Ind,:),Blade.Twist(Ind,:), Plots);

fprintf('Final AEP obtained: %d \n', AEP)
fprintf('Lambda for Chord: %d \n ', (Blade.Optimization.Lambda(Ind)))
fprintf('Pitch offset: %d \n', (Blade.Optimization.Pitch(Ind)))


%Save_Function(Operation,Blade,Ind,Data);



