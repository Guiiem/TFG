function Blade = Creator_Validation (data)

%We will use the linear extrapolation. Set the minimum lambda
R = data.r_tip; %Blade radius
Efficiency = data.cl_array./data.cd_array;
[~, index] = (max(Efficiency)); %Optimal lift coefficient
Cl_opt = mean(data.cl_array(index));
a_opt = deg2rad(mean(data.alpha(index))); %Optimal angle of attack
N_bl = data.numb; %Number of blades
mu = data.rad; %Position of each node


%Lambda range
Lambda = 2.5;

len = length(Lambda);

Blade.Chord = zeros(len,data.no_genes);
Blade.Twist = zeros(len,data.no_genes);


%Generate all chord and twist possible values
for i = 1:len
    for j = 1:data.no_genes
        Blade.Chord(i,j) = R*2*pi*8/(9*Lambda(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*mu(j))^2*(1+2/(9*(Lambda(i)*mu(j))^2))^2)); %Optimal chord distribution
        %Blade.Chord(i,j) =
        %16*pi*R/(9*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*mu(j)+2/(9*Lambda(i)*mu(j)))^2));
        %Small wind turbine book
        %Blade.Chord(i,j) = (R*2*pi*8 / (Lambda(i)^2 *N_bl *Cl_opt *9*0.8)) *(2 - Lambda(i)*mu(j) / (Lambda(i)*0.8)); %Lineal chord distribution
        Blade.Twist(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu(j)*(1+2/(9*(Lambda(i)*mu(j))^2))))-a_opt);
    end
    
    Blade.Chord_Lineal(i,:) = ChordLinealization(Blade.Chord(i,:),mu); %Chord lineal
    Blade.Chord_Lineal_Limited(i,:) = ChordLimitation(Blade.Chord_Lineal(i,:)); %Chord lineal + limited to a certain value
    Blade.Chord_Limited(i,:) = ChordLimitation(Blade.Chord(i,:)); %Chord limited to a certain value. The rest stands as the optimal geometry
end


%We have the values at the center of each blade elements, we can also
%compute the chord and the twist at the boundary of each element (we will
%use it to have the complete geomtery of the blade)

mu_2 = data.rad_boundary/R;

for i = 1:len   
    for j = 1:data.no_genes+1
        Blade.Chord_Boundary(i,j) = R*2*pi*8/(9*Lambda(i)*N_bl*Cl_opt*sqrt(4/9+(Lambda(i)*mu_2(j))^2*(1+2/(9*(Lambda(i)*mu_2(j))^2))^2));
        Blade.Twist_Boundary(i,j) = rad2deg(atan((2/3)/(Lambda(i)*mu_2(j)*(1+2/(9*(Lambda(i)*mu_2(j))^2))))-a_opt);
    end
end

Blade.Radius_Boundary = data.rad_boundary;
Blade.Lambda = Lambda;

%Eliminate the blades that have chord solidity higher than one
Cs = 0.5/pi*N_bl*Blade.Chord(:,1)./(mu(1)*R);
ind = 0;
for i = 1:len
    if Cs(i) >= 1
        ind = i;
    end
end

Blade.Chord_Limited = Blade.Chord_Limited((ind+1:len),:);
Blade.Chord = Blade.Chord((ind+1:len),:);
Blade.Chord_Lineal = Blade.Chord_Lineal((ind+1:len),:);
Blade.Twist = Blade.Twist((ind+1:len),:);
Blade.Chord_Boundary = Blade.Chord_Boundary((ind+1:len),:);
Blade.Twist_Boundary = Blade.Twist_Boundary((ind+1:len),:);
Blade.Lambda = Blade.Lambda(ind+1:len);


        

%cs = 0.5/pi*Nb*Chord./(mu*D/2);





end