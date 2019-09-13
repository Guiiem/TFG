function [CP,CT,Cp] = PowerCoef(Data,Chord,Twist)
%PowerCoef: The power coefficient for a node will be computed for a
%tip speed ratio and a pitch angle given. 

lambda = Data.lambda_power;
tau = Twist;
theta = 0;
mu = Data.rad;
Nb = Data.numb;
cs = 0.5/pi*Numb*Chord./mu; %Local solidity
D = Data.r_tip*2; %Rotor diameter
r = mu*D/2; %Radial positions
Delta_r = r(2)-r(1);
delta = 5.e-3;

%Start values
a_calc = 1/3; %Start considering the optimal axial induction factor 
a = -1/3;
ite_max = 50; %Maximum number of iterations 
CP = 0; %Power coefficient
CT = 0; %Thrust coefficient

torque = 0;

for j=1:length(r) %We will iterate for each node
    
    for i=1:ite_max

        %Induction factors
        a = abs(a-a_calc)*0.5; %Average for the new iteration if the previous value is acceptable
        
        if i==1
            ap_calc = a_calc*(1-a_calc)/(lambda^2 * mu(j)^2); %Computation of the tangencial induction factor
        end
        ap = ap_calc;
            
        %Finding the angle of attack
        phi = atan( (1-a) / (lambda * mu(j) * (1+ap))); %Computation of the angle phi (from rotor plane to relative velocity)
        beta = theta + tau(j); %Resultant of the pitch and the twist
        alpha = phi*(180/pi) - beta; %Angle of attack
        
        %Interpolation of the coefficients of lift and drag for the found alpha
        [Cl,Cd] = Coefficients (alpha);

        %Calculate the induction factors with these values
        A = (cs(j) / (4*sin(phi)^2) )*( Cl*cos(phi) + Cd*sin(phi) );
        a_calc = A / (1+A);

        AP = (cs(j) / (4*sin(phi)*cos(phi)) ) * ( Cl*sin(phi) - Cd*cos(phi) );
        ap_calc = AP / (1-AP);

        %Check the convergency 

        if(abs(a_calc-a) < delta && i>2) 
            i
            break;
        end
    end

    a = a_calc; 
    ap = ap_calc;
    
    %%Nou càlcul de coefficients
%     Ut = (1-a)/sin(phi);
%     
%     delthr = Nb*Ut^2*c(j)*Delta_r/pi; 
%     deltor = delthr*r(j)*(Cl*sin(phi) - Cd*cos(phi));
%     
%     CT = CT + delthr;
%     torque = torque + deltor;
    
    
    
    
    
    
    
    
%%Vell càlcul de coefficients
    Cp(j) = 4*a*(1-a)^2; % Coefficient of power of the node
    
    if (Cp(j)>16/27)
        Cp (j) = 0.005;
    end
    
    if (Cp(j)<0)
        Cp(j) = 0;
    end
        
    CP = CP + Cp(j) * (2*pi*r(j)*Delta_r) / (pi*(D/2)^2); % Total coefficient of power
    
    Ct = 4*a*(1-a); %Coefficient of thrust of the node
    CT = CT + Ct * (2*pi*r(j)*Delta_r) / (pi*(D/2)^2); % Total coefficient of thrust
end

% CP = torque*lambda;

end

