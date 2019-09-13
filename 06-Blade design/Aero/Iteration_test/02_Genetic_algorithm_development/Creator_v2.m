function Blade = Creator_v2 (data)

Blade.Chord = zeros(data.no_indiv,data.no_genes);
Blade.Twist = zeros(data.no_indiv,data.no_genes);

Blade.fitness = zeros(data.no_indiv,2);
Blade.dom = true(data.no_indiv,1);
Blade.age = ones(data.no_indiv,1);

%Select the range of lambda to generate the random values
%Would be a good idea to separate these values for the convergence ones
%Lambda_min = data.lambda_range(1);
%Lambda_max = data.lambda_range(length(data.lambda_range));
Lambda_min = 2.5;
Lambda_max = 10;

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
    Lambda = Lambda_min + (Lambda_max-Lambda_min)*rand;
    for j = 1:data.no_genes
        Blade.Chord(i,j) = R*2*pi*8/(9*Lambda*N_bl*Cl_opt*sqrt(4/9+(Lambda*mu(j))^2*(1+2/(9*(Lambda*mu(j))^2))^2));
        Blade.Twist(i,j) = rad2deg(atan((2/3)/(Lambda*mu(j)*(1+2/(9*(Lambda*mu(j))^2))))-a_opt);
    end
end



end