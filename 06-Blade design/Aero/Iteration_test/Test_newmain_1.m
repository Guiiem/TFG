clc 
clear variables
close all

%%Test main for the new iterative procedure

%%Inputs (all this information should latter be moved to an outer input
%%file
Nb = 3; %Number of blades
D = 1.5; %Rotor diameter
Rho = 1.225;
cmax = 0.25; %Maxchord
cmin = 0.05; %Minchord
N = 15; %Number of blade elements
Umax = 20; %Cut out wind speed
twistmax = 0.1;
twistmin = 0;
rhub = 0.15; 
Delta = 1e-4; %Convergence criteria

%Vectors of radius, chord and twist
Delta_r = (D/2-rhub)/N;
r = (rhub+Delta_r/2):((D/2-rhub)/N):(D/2-Delta_r/2);
Mu = r./(D/2);
c = fliplr(cmin:((cmax-cmin)/(N-1)):cmax);
twist = fliplr(twistmin:((twistmax-twistmin)/(N-1)):twistmax);
Sigma_r = (Nb*c)./(2*pi*r); %Chord solidity

%Curve torque-rpm generator
TGen = [0 25 50 75 100];
OmegaGen = [0 110 200 290 400];

%Rang de lambda a analitzar
Lambda = 0.5:0.1:12;


for i=1:length(Lambda)
        %De l'altre codi
        [CP(i),CT(i),Cp(i,:)] = PowerCoef(0, twist, Lambda(i), Mu, Sigma_r, Delta, r, Delta_r, D, c, Nb);
end

%CP vs lambda
figure('Name','Cp vs Lambda')
plot(Lambda,CP)

%Power vs omega at eaach wind speed
for i=1:Umax
    for j=1:length(Lambda)
        P(i,j) = CP(j)*0.5*Rho*pi*(D/2)^2*i^3; %Aerodynamic power for each (wind speed,lambda)
        Omega(i,j) = Lambda(j)*i/(D/2); %Rotor speed for each (wind speed, lambda)
        Q(i,j) = P(i,j)/Omega(i,j); %Aerodynamic torque for each (wind speed, lambda)
    end
end

figure('Name','Torque vs Omega for each lambda')
for i=1:length(Lambda)
    plot(Omega(:,i),Q(:,i));
    hold on
end

figure('Name','Torque vs Omega for each wind speed')
for i=1:Umax
    plot(Omega(i,:),Q(i,:));
    hold on
    legendInfo{i}=['U: ' num2str(i)];
end
hold on
plot(OmegaGen,TGen)
legend(legendInfo,'Generator curve');



Q(isnan(Q))=0;

%On es creuen a u=7, per exempe?
for i=1:Umax
    [Qeq(i) Omegaeq(i)] = Intersect(Q(i,:),Omega(i,:),TGen,OmegaGen,i);
end

figure
plot(Omegaeq,Qeq)




