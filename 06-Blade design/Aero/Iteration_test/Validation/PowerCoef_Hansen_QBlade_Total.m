function Power_Curve = PowerCoef_Hansen_QBlade_Total(Options, Power_Curve,Data,Chord_org,Twist_org,m)

lambda = Data.lambda_power;
mu = Data.rad;
Nb = Data.numb;
D = Data.r_tip*2; %Rotor diameter
R_hub = Data.r_hub; %Hub radius
%cs = 0.5/pi*Nb*Chord./(mu*D/2); %Local solidity
r_org = Data.rad_boundary; %Radial positions
Delta_r_org = r_org(2)-r_org(1);
delta = 0.001;

%Lift and drag coefficients
aoa_int=Data.alpha; % Alpha array for Cl and Cd data
aoa_max=max(aoa_int);% Max. angle of attack
aoa_min=min(aoa_int); % Min. angle of attack
i_max=length(aoa_int); % Length of alpha array
delalpha=Data.del_alpha; % Increment in aoa
cl_array=Data.cl_array(:,1); % Cl array
cd_array=Data.cd_array(:,1); % Cd array
Re_array=Data.re_array; % Re array
num_re=length(Re_array); % length of Re array
Re_log=Data.re_log; % log(Re) array

%New inputs
V0 = Data.u_power; %Incoming wind speed [m/s]
omega = lambda*V0/(D/2); %Rotor speed [rad/s] 
visc = 1.647e-5; % Kinematic viscosity of air (m^2/s) needed for Reynolds no.
rho = 1.225; %Air density
rel_fac = 0.5; %Relaxation factor

%Implementation of the sinusoidal spacing:

sections = length(r_org)-1;
deltaangle = pi/(sections+1);
len = D/2 - R_hub;
tot = 0;

for i=1:sections
    ang = i*deltaangle;
    dist(i) = sin(ang);
    tot = tot + sin(ang);
end
thisdelta = len/tot;

for i=1:sections
    Delta_r(i) = thisdelta*dist(i);
end

position = R_hub;

%Start the discretization
for i=1:sections
    
    %position marks the center of each element
    position = position + Delta_r(i)/2;
    r(i) = position;
    
    %Interpolate the value of chord and twist of the element
    if position==r_org(sections)
        Chord(i) = Chord_org(sections);
        Twist(i) = Twist_org(sections);
    else
        i_rad = find(r_org>position,1);
        factor = (position-r_org(i_rad-1))/Delta_r_org;
        Chord(i) = Chord_org(i_rad-1)*(1-factor)+Chord_org(i_rad)*factor;
        Twist(i) = Twist_org(i_rad-1)*(1-factor)+Twist_org(i_rad)*factor;
    end
    position = position + Delta_r(i)/2;

end

tau = Twist;
mu = r./(D/2);    
    






for i=1:sections
    %Initialize induction factors
    a = 0; a_calc = 0; a_prev = 0;
    ap = 0; ap_calc = 0;
    for j=1:300
        
        %Relaxation factor
        if j>100
            a_new = a_calc*1/4 + a*1/2 + a_prev*1/4;
        elseif j>10 
            a_new = a_calc*rel_fac + (1-rel_fac)*a;
        else
            a_new = a_calc;
        end
        
        ap_new =ap_calc;
        
        a_prev = a; ap_prev = ap;
        a = a_new; ap = ap_new;
        
        %Compute flow angle
        phi = atan( (1-a)*V0/((1+ap)*omega*r(i)) );
        aoa = phi*180/pi - tau(i);
        
        cosphi = cos(phi);
        sinphi = sin(phi);
        
        %Find aerodynamic coefficients
        Vt = V0*(1-a)/sinphi; %Total or relative velocity
        Rei = Vt*Chord(i)*rho/visc;
        
        if Rei <= Re_array(1) % For Re less than minimum measured
              cl_int=cl_array; %(:,1);
              cd_int=cd_array; %(:,1);
        elseif Rei >= Re_array(num_re) % For Re greater than max. measured
              cl_int=cl_array; %(:,num_re);
              cd_int=cd_array; %(:,num_re);
        else % Re is in the tabulated range
            Re_lg=log(Rei);
            j_Re=find(Re_log>Re_lg,1);
            factor_Re = (Re_lg - Re_log(j_Re-1))/(Re_log(j_Re) - Re_log(j_Re-1));
            %Now interpolate in log10(Re).
            cl_int=cl_array; %(:,j_Re-1)*(1-factor_Re)+cl_array(:,j_Re)*factor_Re;
            cd_int=cd_array; %(:,j_Re-1)*(1-factor_Re)+cd_array(:,j_Re)*factor_Re;
        end
        if (aoa >= aoa_max) % Use generic high-alpha equations
            Cl=sind(2*aoa);
            Cd=(sind(aoa))^2;
        elseif (aoa <= aoa_min) % Penalise blades for low alpha
            C_a=0;
            C_adash=0;
            break
        else
            i_up=find(aoa_int>aoa,1);
            factor= (aoa-aoa_int(i_up-1))/delalpha;
            Cl=cl_int(i_up-1)*(1-factor)+ cl_int(i_up)*factor;
            Cd=cd_int(i_up-1)*(1-factor)+ cd_int(i_up)*factor;
        end
        
        Cn = Cl*cosphi + Cd*sinphi;
        Ct = Cl*sinphi - Cd*cosphi;
        
        %Computation of solidity
        cs = Chord(i)*abs(cos(deg2rad(Twist(i))))*Nb/(2*pi*r(i)); %Modified from QBlade
        
        %Prandtl's tip loss correction
        f = sinphi;
        g = (D/2 - r(i))/r(i);
        F = 2/pi * acos(exp(-Nb/2 * abs(g/f)));
                
        %Converge of the first iterations
        if ~isreal(F) || Options.TipLoss == false
            F=1;
        end
        
        %Root loss
        if Options.RootLoss == 1
            g = (r(i)-R_hub)/r(i);
            Froot = 2/pi * acos(exp(-Nb/2 * abs(g/f)));
            F = F*Froot;
        end
        
        %Computation of the local trhust coefficient
        CT = cs*(1-a)^2 * Cn /(sinphi^2);
        if CT <= 0.96*F
            a_calc = 1/((4*F*sinphi^2)/(cs*Cn) +1);
        else 
            a_calc = (18*F-20-3*sqrt(abs(CT*(50-36*F)+12*F*(3*F-4))))/(36*F-50);
        end
      
        ap_calc = 0.5*(sqrt(abs(1+(4*a_calc*(1-a_calc)/((lambda*mu(i))^2))))-1);
        diff = min(abs(a_calc-a),abs(ap_calc-ap));
        
        if (j>2 && abs(a_calc-a)<delta && abs(ap_calc-ap)<delta)
            a = a_calc;
            ap = ap_calc;
            break
        end
    end
    
    %Normal and tangencial dimensional components of the blade section load
    Vt = sqrt((V0*(1-a_calc))^2+(V0*lambda*mu(i)*(1+ap_calc))^2);
    Pn(i) = Cn*0.5*rho*Vt^2*Chord(i);
    Pt(i) = Ct*0.5*rho*Vt^2*Chord(i);
    
    
    
    Power_Curve.a(m,i) = a;
    Power_Curve.ap(m,i) = ap;
    Power_Curve.NumItrt(m,i) = j;
    Power_Curve.aoa(m,i) = aoa;
    Power_Curve.Cl(m,i) = Cl;
    Power_Curve.Cd(m,i) = Cd;
    Power_Curve.Re(m,i) = Rei;     
    Power_Curve.Vtot(m,i) = Vt;
    Power_Curve.Pn(m,i) = Pn(i);
    Power_Curve.Pt(m,i) = Pt(i);
    Power_Curve.Phi(m,i) = rad2deg(phi);
    Power_Curve.F(m,i) = F;
    Power_Curve.Cn(m,i) = Cn;
    Power_Curve.Ct(m,i) = Ct;
    Power_Curve.diff(m,i) = diff;

     
    
          
    
end

%Integrate the total shaft torque
M = 0;
for i=1:length(r)-1
    A = (Pt(i+1)-Pt(i))/(r(i+1)-r(i));
    B = (Pt(i)*r(i+1) - Pt(i+1)*r(i))/(r(i+1)-r(i));
    
    Mi(i) = 1/3 * A*(r(i+1)^3 -r(i)^3) + 1/2*B*(r(i+1)^2-r(i)^2);
    Power_Curve.TorqueContribution(m,i) = Mi(i);   
    
    M = Nb*Mi(i) + M;
end

Power = M*omega;

Power_Curve.Cp(m) = Power/(0.5*rho*V0^3*(D/2)^2*pi);
Power_Curve.Power(m) = Power;
Power_Curve.Torque(m) = M;
Power_Curve.Omega(m) = omega*60/(2*pi);


%Calculation of the power coefficient
% power = 0;
% for i=1:sections
%     power = power + r(i)*Power_Curve.Pt(m,i)*Delta_r(i)/V0^2;
% end
% 
% Power_Curve.Power(m) = power*Nb*lambda/(D/2);
% windenergy = 0.5*pi*(D/2)^2*rho;
% Power_Curve.Cp(m) = power/windenergy;


Power_Curve.rad = r;
Power_Curve.lambda(m) = lambda;