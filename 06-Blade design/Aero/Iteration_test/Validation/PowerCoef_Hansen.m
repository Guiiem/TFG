function Power_Curve = PowerCoef_Hansen(Power_Curve,Data,Chord,Twist,m)

lambda = Data.lambda_power;
tau = Twist;
mu = Data.rad;
Nb = Data.numb;
D = Data.r_tip*2; %Rotor diameter
cs = 0.5/pi*Nb*Chord./(mu*D/2); %Local solidity
r = mu*D/2; %Radial positions
Delta_r = r(2)-r(1);
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
rel_fac = 0.7; %Relaxation factor

for i=1:length(r)
    %Initialize induction factors
    a = 0; a_calc = 0;
    ap = 0; ap_calc = 0;
    for j=1:1000
        a = a*rel_fac + a_calc*(1-rel_fac);
        ap = ap*rel_fac + ap_calc*(1-rel_fac);
        
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
        
        %Glauert's correction
        F = 1;
        K = 4*F*sinphi^2/(cs(i)*Cn);
        ac = 0.2;
        
        a_calc = 1/( (4*sinphi^2/(cs(i)*Cn)) +1);
        if a_calc > ac && (K*(1-2*ac)+2)^2 + 4*(K*ac^2-1)>=0
            a_calc = 0.5*(2+K*(1-2*ac)-sqrt( (K*(1-2*ac)+2)^2 + 4*(K*ac^2-1)));
        end
        
        ap_calc = 1/( (4*sinphi*cosphi/(cs(i)*Ct)) -1);
        
        if (j>2 && abs(a_calc-a)<delta && abs(ap_calc-ap)<delta)
            a = a_calc;
            ap = ap_calc;
            break
        end
    end
    
    %Normal and tangencial dimensional components of the blade section load
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
Power_Curve.Omega(m) = omega;




Power_Curve.lambda(m) = lambda;