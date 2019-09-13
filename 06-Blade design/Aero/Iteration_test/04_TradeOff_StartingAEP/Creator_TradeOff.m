function Blade = Creator_TradeOff (data)

Blade.Chord = zeros(data.no_indiv,data.no_genes);
Blade.Twist = zeros(data.no_indiv,data.no_genes);

Blade.fitness = zeros(data.no_indiv,2);

%Select the range of lambda to generate the random values
%Would be a good idea to separate these values for the convergence ones
%Lambda_min = data.lambda_range(1);
%Lambda_max = data.lambda_range(length(data.lambda_range));
Lambda_min = 1.5;
Lambda_max = 4;

Lambda = Lambda_min:((Lambda_max-Lambda_min)/data.no_indiv):Lambda_max;

%Obtain the necessary parameters to compute the optimal chord and twist
R = data.r_tip; %Blade radius
Efficiency = data.cl_array./data.cd_array;
[~, index] = (max(Efficiency)); %Optimal lift coefficient
Cl_opt = mean(data.cl_array(index));
a_opt = deg2rad(mean(data.alpha(index))); %Optimal angle of attack
N_bl = data.numb; %Number of blades
mu = data.rad; %Position of each node


%Minimum and maximum chord and twist are no longer used for generation

for i = 1:data.no_indiv    % Generate the first population
    for j = 1:data.no_genes
        Blade.Chord(i,j) = R*2*pi*8/(9*Lambda(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*mu(j))^2*(1+2/(9*(Lambda(i)*mu(j))^2))^2));
        Blade.Twist(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu(j)*(1+2/(9*(Lambda(i)*mu(j))^2))))-a_opt);
    end
end

%We have the values at the center of each blade elements, we can also
%compute the chord and the twist at the boundary of each element (we will
%use it to have the complete geomtery of the blade)

mu_2 = data.rad_boundary/R;

for i = 1:data.no_indiv    % Generate the first population
    for j = 1:data.no_genes+1
        Blade.Chord_Boundary(i,j) = R*2*pi*8/(9*Lambda(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*mu_2(j))^2*(1+2/(9*(Lambda(i)*mu_2(j))^2))^2));
        Blade.Twist_Boundary(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu_2(j)*(1+2/(9*(Lambda(i)*mu_2(j))^2))))-a_opt);
    end
end

Blade.Radius_Boundary = data.rad_boundary;




end