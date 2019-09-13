%%Script for testing BEMT functions
clear variables
close all

addpath ..\01_Non_iterative_development\;

Data = blade_opt_setup;

%Optimize for lambda = 6;
R = Data.r_tip; %Blade radius
Efficiency = Data.cl_array./Data.cd_array;
[~, index] = (max(Efficiency)); %Optimal lift coefficient
Cl_opt = mean(Data.cl_array(index));
a_opt = deg2rad(mean(Data.alpha(index))); %Optimal angle of attack
N_bl = Data.numb; %Number of blades
mu = Data.rad; %Position of each node

Lambda = 6;

Chord = zeros(Data.no_genes,1);
Twist = zeros(Data.no_genes,1);

for i=1:Data.no_genes
    Chord(i) = R*2*pi*8/(9*Lambda*N_bl*Cl_opt*sqrt(4/9+(Lambda*mu(i))^2*(1+2/(9*(Lambda*mu(i))^2))^2));
    Twist(i) = rad2deg(atan((2/3)/(Lambda*mu(i)*(1+2/(9*(Lambda*mu(i))^2))))-a_opt);
end

%Compute with book's function
Data.lambda_power = Lambda;
[Cp1, Cq1] = power_calc_opt(Data,Chord,Twist,0);

%%Compute with new function
[Cp2, Cq2] = PowerCoef_v2(Data,Chord,Twist);




