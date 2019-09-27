function [cp, cq, ct] =power_calc_opt(data,chord,twist, iprint)
% Extract information from data structure
aoa_int=data.alpha; % Alpha array for Cl and Cd data
aoa_max=max(aoa_int);% Max. angle of attack
aoa_min=min(aoa_int); % Min. angle of attack
i_max=length(aoa_int); % Length of alpha array
delalpha=data.del_alpha; % Increment in aoa
r_tip=data.r_tip; % Blade tip radius
cl_array=data.cl_array; % Cl array
cd_array=data.cd_array; % Cd array
Re_array=data.re_array; % Re array
num_re=length(Re_array); % length of Re array
Re_log=data.re_log; % log(Re) array
Numb=data.numb; % Number of blades
U0=data.u_power; % Wind speed of rpower calculation
lambda = data.lambda_power; % TSR for power calculations
rad=data.rad;  % Radius of blade elements
nbes = data.no_genes; % No. blade elements

visc = 1.5e-5; % Kinematic viscosity of air (m^2/s) needed for Reynolds no.
tol = 5.e-3;   % Convergence tolerance for BE analysis
delr=rad(2)-rad(1); % Determine width of blade elements
lamr = lambda*rad; % Local speed ratio
smallf = 0.5*Numb*(1./rad-1);
thrust = 0.0; torque = 0.0;
a = 0.3; % Initialise a
sigma = 0.5/pi*Numb*chord./rad;  % Local solidity
a_c = 1/3.0; % Parameter for high thrust correction
Re = U0*chord.*r_tip./visc ; % The Reynolds number

heading_screen=...
    '   Radius  iter.  aoa       a       Cl       Cd       deltor     Re \n';
format_screen='  %7.4f %3d  %7.2f  %7.3f  %7.3f  %8.5f  %8.5f  %8.3e \n';
if (iprint == 1); fprintf(heading_screen); end;

for i = 1: nbes
      adash=0.0;
      Ut=sqrt((1-a)^2+lamr(i)^2);
      for j = 1:20 % Limit the  number of iterations per BE to 20
              phi = atan((1 - a)/(1 + adash)/lamr(i));
              cosphi = cos(phi);
              sinphi = sin(phi);
              bigF = 2*acos(exp(-smallf(i)/sinphi))/pi; % Equation 5.1 -> Prandtl's tip loss factor
%             bigF = 1.0; %(max value)
              aoa = phi*180.0/pi - twist(i); % Find aoa
              Ut=(1-a)/sinphi;
              Rei = Ut*Re(i);  % The Reynolds number
              
              %Interpolate values of lift and drag coefficients
              if Rei <= Re_array(1) % For Re less than minimum measured
                    cl_int=cl_array(:,1);
                    cd_int=cd_array(:,1);
              elseif Rei >= Re_array(num_re) % For Re greater than max. measured
                    cl_int=cl_array(:,num_re);
                    cd_int=cd_array(:,num_re);
              else % Re is in the tabulated range
                    Re_lg=log(Rei);
                    j_Re=find(Re_log>Re_lg,1);
                    factor_Re = (Re_lg - Re_log(j_Re-1))/(Re_log(j_Re) - Re_log(j_Re-1));
                    % Now interpolate in log10(Re).
                    cl_int=cl_array(:,j_Re-1)*(1-factor_Re)+cl_array(:,j_Re)*factor_Re;
                    cd_int=cd_array(:,j_Re-1)*(1-factor_Re)+cd_array(:,j_Re)*factor_Re;
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
              
              
              C_a = Cl*cosphi + Cd*sinphi; %Resultant in the circumferencial direction. Parameter defined by us (eq. 3.10)
              Y1 = 4*bigF*sinphi^2/(sigma(i)*C_a); %Intermediate functions for the axial induction factor calculation
              
              if (Y1 >= 2.0 ) % In the low thrust region
 %             if (Y1 >= 1.0 ) % In the low thrust region
                    newa = 0.5*(2 + Y1 - sqrt(Y1*(4*(1-bigF) + Y1)))/(1+bigF*Y1);
              else % Use Equation (2.20) for the high thrust region
                    c_t = 1 -2*a_c*bigF;
                    newa = 1.0 + 0.5*(Y1*c_t-sqrt((Y1*c_t+2)^2 -4*(1-a_c^2*bigF*Y1)));
              end
              diffa = abs(a - newa); 
              
              %Provisional method to prevent problems in the first iteration 
              if newa>=1
                  newa = 0.2;
              end
              
              a = real(newa);
              C_adash = Cl*sinphi - Cd*cosphi; % Resultant in circumferential direction (Parameter defined by us (eq. 3.11))
              Y2 = 4*bigF*sinphi*cosphi/sigma(i)/C_adash;
              adash = real((1-a)/((1-a*bigF)*Y2 +a-1));
              if j > 2 && diffa < tol*a; break;  end
      end
      delthr = Numb*Ut*Ut*chord(i)*delr/pi;  % Equn (3.10)
      deltor = delthr*rad(i)*C_adash;   % Equn (3.11)
      delthr = delthr*C_a; % Completing equn (3.10)
      thrust = thrust + delthr;
      torque = torque + deltor; % Sum the rotor torque
      if (iprint == 1)
             fprintf(format_screen, rad(i), j, aoa, a, Cl, Cd, deltor, Rei)
      end
cp(i) = torque*lambda;
if cp > 0.593; cp=0.005;end

cq(i) = torque;
ct(i) = thrust;

%Informació a guardar: Cp, NumItrt, Y1, F, f, AoA, Cl, Cd, a, ap, Blade
%solidity, meancl, meanaoa, lambda, cp, ct, cq
end

