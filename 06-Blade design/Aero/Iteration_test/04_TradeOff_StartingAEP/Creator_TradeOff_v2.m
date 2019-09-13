function Blade = Creator_TradeOff_v2 (data)

%V2: chord and twist are not optimized for the same lambda
%Chord is now limited to 20cm

Blade.Chord = zeros(data.no_indiv,data.no_genes);
Blade.Twist = zeros(data.no_indiv,data.no_genes);

Blade.fitness = zeros(data.no_indiv,2);

indiv_para = data.no_indiv^(1/3); %Number of individuals for each parameter

%Select the range of lambda to generate the random values
%Would be a good idea to separate these values for the convergence ones
%Lambda_min = data.lambda_range(1);
%Lambda_max = data.lambda_range(length(data.lambda_range));
Lambda_min = 0.2;
Lambda_max = 5;

Pitch_min = -2;
Pitch_max = 2;


Lambda = Lambda_min:((Lambda_max-Lambda_min)/(indiv_para-1)):Lambda_max;
Pitch = Pitch_min:((Pitch_max-Pitch_min)/(indiv_para)):Pitch_max;

%Obtain the necessary parameters to compute the optimal chord and twist
R = data.r_tip; %Blade radius
Efficiency = data.cl_array./data.cd_array;
[~, index] = (max(Efficiency)); %Optimal lift coefficient
Cl_opt = mean(data.cl_array(index));
a_opt = deg2rad(mean(data.alpha(index))); %Optimal angle of attack
N_bl = data.numb; %Number of blades
mu = data.rad; %Position of each node

%Find the lambda at which the chord is 20 cm
for i=1:length(Lambda)
    Chord = R*2*pi*8/(9*Lambda(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*min(mu))^2*(1+2/(9*(Lambda(i)*min(mu))^2))^2));
    if Chord < 0.2
        Lambda_min_chord = Lambda(i);
        break;
    end
end

%Redefine lambda vector for chord
Lambda_c = Lambda_min_chord:((Lambda_max-Lambda_min_chord)/(indiv_para-1)):Lambda_max;

indiv_para = int16(indiv_para);

%Generate all chord and twist possible values
for i = 1:indiv_para 
    for j = 1:data.no_genes
        Chord(i,j) = R*2*pi*8/(9*Lambda_c(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda_c(i)*mu(j))^2*(1+2/(9*(Lambda_c(i)*mu(j))^2))^2));
        Twist(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu(j)*(1+2/(9*(Lambda(i)*mu(j))^2))))-a_opt);
    end
end

%Create the distribution
for i = 1:indiv_para
    for j = 1:indiv_para
        for k = 1:indiv_para
            ind = i*indiv_para^2+j*indiv_para-indiv_para^2-indiv_para+k;
            Blade.Chord(ind,:) = Chord(i,:);
            Blade.Twist(ind,:) = Twist(j,:)+Pitch(k);
            
            Blade.Optimization.Chord_Lambda(ind) = Lambda_c(i);
            Blade.Optimization.Twist_Lambda(ind) = Lambda(j);
            Blade.Optimization.Pitch(ind) = Pitch(k);
        end
    end 
end
        

%We have the values at the center of each blade elements, we can also
%compute the chord and the twist at the boundary of each element (we will
%use it to have the complete geomtery of the blade)

mu_2 = data.rad_boundary/R;

for i = 1:indiv_para    
    for j = 1:data.no_genes+1
        Chord_Boundary(i,j) = R*2*pi*8/(9*Lambda_c(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda_c(i)*mu_2(j))^2*(1+2/(9*(Lambda_c(i)*mu_2(j))^2))^2));
        Twist_Boundary(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu_2(j)*(1+2/(9*(Lambda(i)*mu_2(j))^2))))-a_opt);
    end
end

%Create the distribution
for i = 1:indiv_para
    for j = 1:indiv_para
        for k = 1:indiv_para
            ind = i*indiv_para^2+j*indiv_para-indiv_para^2-indiv_para+k;
            Blade.Chord_Boundary(ind,:) = Chord_Boundary(i,:);
            Blade.Twist_Boundary(ind,:) = Twist_Boundary(j,:)+Pitch(k);
        end
    end 
end

Blade.Radius_Boundary = data.rad_boundary;



end