clc
close all
clear variables

addpath ..\01_Non_iterative_development\;
addpath ..\03_PowerCalcTest\;

tic

%%Generation of the Data file to be used in the PowerCoef Function
Data = blade_opt_setup;

%%Select which plots you want to be printed or not
Plots = PlotSelection;

%%Generation of the initial population
Blade = Creator_v2(Data);

%Compute the AEP for each individual of the initial generation
for i=1:Data.no_indiv
    
    [Blade.AEP(i), Data, ~] = BigSolver (Data, Blade.Chord(i,:), Blade.Twist(i,:), Plots);
    
end

%%Start the optimization process through all the generations 
for gen=1:Data.no_gen  
        
    %Find the fittest members of the current population
    [i_max_fit, max_fit, Blade, Data] = Fittest_v2(gen, Data, Blade);
    Data.max_fit(gen,:) = max_fit;
    
    %Breed population, except for the last generation
    if gen < Data.no_gen
        [Blade, Data] = Breed_v2 (gen, i_max_fit, max_fit, Data, Blade, Plots);
    end
    
end

%Compute the operation curves for the "best blade"
Plots = PlotSelectionFinal; %Now we want to see the final plots

[~, ~, Operation] = BigSolver(Data, Blade.Chord(i_max_fit,:),Blade.Twist(i_max_fit,:), Plots);




simulation_time = toc;

