%%First approach to the needed beam

Max_stress = 1200*10^6; %[Pa]

%Geometric data
a = [4 8 10 10]*10^-3; %[m]
b = [2 7 8 4]*10^-3; %[m]

%Moment of area
for i=1:3
    Iz(i) = a(i)^4/12-b(i)^4/12;
    y(i) = a(i)/2;
end

Iz(4) = a(4)^4/12-pi*b(4)^4/4;
y(4) = a(4)/2;

%Maximum load calculation
Mz = Max_stress.*Iz./y;

leg = categorical({'S1: a=4, b=2', 'S2: a=8, b=7', 'S3: a=10, b=8', 'S4: a=10, r=4'});
bar(leg,Mz)
ylabel('Maximum bending moment [Nm]')

