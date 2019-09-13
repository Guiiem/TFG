function Data = Viterna_Extrapolation(Data)

Cl = Data.cl_array;
Cd = Data.cd_array;
Alpha = Data.alpha;
Polar_len = length(Alpha);
Num_Re = Data.num_re;

AR = 0.55/0.12; %Approximated Aspect Ratio
Cd_max = 1.11 + 0.18*AR;

%First constant calculation
A1 = Cd_max/2;
B1 = Cd_max;

for i = 1:Num_Re
    %Parameters for compute the secondary constants
    Cl_s(i) = Cl(Polar_len,i);
    Cd_s(i) = Cd(Polar_len,i);
    AoA_s = deg2rad(Alpha(Polar_len));
    
    %Secondary constants calculation
    A2(i) = (Cl_s(i)-Cd_max*sin(AoA_s)*cos(AoA_s))*(sin(AoA_s)/(cos(AoA_s)^2));
    B2(i) = Cd_s(i) - (Cd_max*sin(AoA_s)^2)/cos(AoA_s);
    
    %Extension of alpha region 
    AoA_s = round(rad2deg(AoA_s),0);
    AoA_range = 90-AoA_s;
    
    %Create variables
   % Cl_new = zeros(AoA_range*2+Polar_len,Num_Re);
   % Cd_new = zeros(AoA_range*2+Polar_len,Num_Re);

    %First, create the negative curve
    for j=1:AoA_range
        Alpha_new(j) = -91+j;    
        Cl_new(j,i) = A1*sind(2*Alpha_new(j))+A2(i)*((cosd(Alpha_new(j))^2)/sind(Alpha_new(j)));
        Cd_new(j,i) = B1*sind(Alpha_new(j))^2 + B2(i)*cosd(Alpha_new(j));             
    end
    
    %Then, add the polar region we already had
    for j=1:Polar_len
        Alpha_new(j+AoA_range) = Alpha(j);
        Cl_new(j+AoA_range,i) = Cl(j,i);
        Cd_new(j+AoA_range,i) = Cd(j,i);   
    end
    
    %Finally, add the positive region
    for j=1:AoA_range
        Alpha_new(j+AoA_range+Polar_len) = AoA_s+j;
        Cl_new(j+AoA_range+Polar_len,i) = A1*sind(2*Alpha_new(j+AoA_range+Polar_len))+A2(i)*((cosd(Alpha_new(j+AoA_range+Polar_len))^2)/sind(Alpha_new(j+AoA_range+Polar_len)));
        Cd_new(j+AoA_range+Polar_len,i) = B1*sind(Alpha_new(j+AoA_range+Polar_len))^2 + B2(i)*cosd(Alpha_new(j+AoA_range+Polar_len));                
    end
    
end
%Extract the new polars
Data.cl_array = Cl_new;
Data.cd_array = Cd_new;
Data.alpha = Alpha_new;   



%     plot(Alpha_new,Cl_new(:,i))
%     hold on
%     plot(Alpha_new,Cd_new(:,i))
    
    

%Extrat extrapolation data for further extrapolations
Data.Viterna.Cl_s = Cl_s;
Data.Viterna.Cd_s = Cd_s;
Data.Viterna.A1 = A1;
Data.Viterna.A2 = A2;
Data.Viterna.B1 = B1;
Data.Viterna.B2 = A2;


end 