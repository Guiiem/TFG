function Power_Curve = PowerCoef_v2_WithTipLosses(Power_Curve,Data,Chord,Twist,m)
%PowerCoef: The power coefficient for a node will be computed for a
%tip speed ratio and a pitch angle given. 

lambda = Data.lambda_power;
tau = Twist;
theta = 0;
mu = Data.rad;
Nb = Data.numb;
D = Data.r_tip*2; %Rotor diameter
cs = 0.5/pi*Nb*Chord./(mu*D/2); %Local solidity
r = mu*D/2; %Radial positions
Delta_r = r(2)-r(1);
delta = 0.002;

%Lift and drag coefficients
aoa_int=Data.alpha; % Alpha array for Cl and Cd data
aoa_max=max(aoa_int);% Max. angle of attack
aoa_min=min(aoa_int); % Min. angle of attack
i_max=length(aoa_int); % Length of alpha array
delalpha=Data.del_alpha; % Increment in aoa
cl_array=Data.cl_array(:,3); % Cl array
cd_array=Data.cd_array(:,3); % Cd array
%cd_array = zeros(length(cl_array),1); %Drag neglection [TEST]
Re_array=Data.re_array; % Re array
num_re=length(Re_array); % length of Re array
Re_log=Data.re_log; % log(Re) array


%Start values
a_calc = 1/3; %Start considering the optimal axial induction factor 
a = 1/3;
ite_max = 5000; %Maximum number of iterations 
CP = 0; %Power coefficient
CT = 0; %Thrust coefficient

torque = 0;

for j=1:length(r) %We will iterate for each node
    
    a_calc = 1/3; %Start considering the optimal axial induction factor 
    a = 1/3;
    
    for i=1:ite_max
        
        

        %Induction factors
        a = a*0.5+a_calc*0.5; %Average for the new iteration if the previous value is acceptable
        %a = a_calc;
        
        if i==1   %We will add the ideal phi as well too 
            ap_calc = a_calc*(1-a_calc)/(lambda^2 * mu(j)^2); %Computation of the tangencial induction factor
        end
        ap = ap_calc;
            
        %Finding the angle of attack
        phi = atan( (1-a) / (lambda * mu(j) * (1+ap))); %Computation of the angle phi (from rotor plane to relative velocity)
        beta = theta + tau(j); %Resultant of the pitch and the twist
        aoa = phi*(180/pi) - beta; %Angle of attack
        
         %Interpolate values of lift and drag coefficients
            
              if (aoa >= aoa_max) % Use generic high-alpha equations
                    Cl=sind(2*aoa);
                    Cd=(sind(aoa))^2;
              elseif (aoa <= aoa_min) % Penalise blades for low alpha
                    Cl=-sind(2*aoa);
                    Cd=-((sind(aoa))^2);
              else
                     i_up=find(aoa_int>aoa,1);
                     factor= (aoa-aoa_int(i_up-1))/delalpha;            
                     Cl=cl_array(i_up-1)*(1-factor)+ cl_array(i_up)*factor;
                     Cd=cd_array(i_up-1)*(1-factor)+ cd_array(i_up)*factor;
              end

        %Calculate the induction factors with these values
        Cx = Cl*cos(phi) + Cd*sin(phi);
        Cy = Cl*sin(phi) - Cd*cos(phi);
        
        %Introduction of the Prandtl's tip loss factor (F=a/ab) (a is the
        %average induction factor and ab the ones at the blades)
        f = Nb*(D/2-r(j))/(2*r(j)*sin(phi));  
        
        if exp(-f)<=-1 || exp(-f)>=1
            F = NaN; %Indicate that the correction could not be applied             
            A = (cs(j) / (4*sin(phi)^2) )*(Cx);
            a_calc = A / (1+A);
            
            AP = (cs(j) / (4*sin(phi)*cos(phi)) ) * ( Cy );
            ap_calc = AP / (1-AP);            
        else        
            F = 2*acos(exp(-f))/pi;            
            Y1 = 4*F*sin(phi)^2./(cs(j)*Cx);            
            %Use Glauert correction for the hight thrust if necessary
            if Y1>= 2
                if(4*Y1*(1-F)+Y1^2>=0)
                    a_calc = (2+Y1-sqrt(4*Y1*(1-F)+Y1^2))/(2*(1+F *Y1));
                    %ap_calc = 1/(((1-a_calc*F)*Y2)/(1-a_calc) -1);
                    Y2 = 4*F*sin(phi)*cos(phi)/(cs(j)*Cy);                    
                    ap_calc = 1/(((1-a_calc*F)*Y2)/(1-a_calc) -1);
                else
                    F = NaN; %Indicate that the correction could not be applied                    
                    A = (cs(j) / (4*sin(phi)^2) )*(Cx);
                    a_calc = A / (1+A);
                    
                    AP = (cs(j) / (4*sin(phi)*cos(phi)) ) * ( Cy );
                    ap_calc = AP / (1-AP);
                end               
            else
                C_t = 1-2/3*F;
                if (Y1*C_t+2)^2-4*(1-(1/3)^2*F*Y1)<0
                    F = NaN;
                    A = (cs(j) / (4*sin(phi)^2) )*(Cx);
                    a_calc = A / (1+A);
                    
                    AP = (cs(j) / (4*sin(phi)*cos(phi)) ) * ( Cy );
                    ap_calc = AP / (1-AP);
                else
                    a_calc = (2+Y1*C_t-sqrt((Y1*C_t+2)^2-4*(1-(1/3)^2*F*Y1)))/2;
                    Y2 = 4*F*sin(phi)*cos(phi)/(cs(j)*Cy);
                    ap_calc = 1/(((1-a_calc*F)*Y2)/(1-a_calc) -1);
                end               

            end

        end
              %Check the convergency 

        if(abs(a_calc-a) < delta && i>2) 
            i;
            break;
        end
    end

    a = a_calc; 
    ap = ap_calc;
    
    Cp(j) = 4*a*(1-a)^2; % Coefficient of power of the node
    
    if (Cp(j)>16/27)
        Cp (j) = 0.005;
    end
    
    if (Cp(j)<0)
        Cp(j) = 0;
    end
    Power_Curve.CpSection(m,j) = Cp(j);
    Power_Curve.NumIterations(m,j) = i;
    Power_Curve.Y1(m,j) = Y1;
    Power_Curve.F(m,j) = F;
    Power_Curve.AoA(m,j) = aoa;
    Power_Curve.Cl(m,j) = Cl;
    Power_Curve.Cd(m,j) = Cd;
    Power_Curve.a(m,j) = a;
    Power_Curve.ap(m,j) = ap;
    CP = CP + Cp(j) * (2*pi*r(j)*Delta_r) / (pi*(D/2)^2); % Total coefficient of power
    Ct = 4*a*(1-a); %Coefficient of thrust of the node
    CT = CT + Ct * (2*pi*r(j)*Delta_r) / (pi*(D/2)^2); % Total coefficient of thrust
end

Power_Curve.BladeSolidity(m,:) = cs;
Power_Curve.MeanCl(m) = mean(Power_Curve.Cl(m,:));
Power_Curve.MeanAoA(m) = mean(Power_Curve.AoA(m,:));
Power_Curve.Lambda(m) = lambda;
Power_Curve.CP(m) = CP;
Power_Curve.CT(m) = CT;
Power_Curve.CQ(m) = CP/lambda;
end

