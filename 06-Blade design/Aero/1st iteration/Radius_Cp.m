clc
clear variables
close all

%This script aims to compute the R and Cp required to match rated power

P = 500; %Rated power
w = 450:50:550; %Rated rotor speed
eta = 0.74; %Generator efficiency

U = 11; %Wind speed at rated power

rho = 1.225; %Air density

R = 0.4:0.05:1.2; %Blade radius
S = pi*R.^2; %Sweept area

Cp = P./(eta*0.5*rho*U^3*S); %Power coefficient at rated power
for i=1:length(w)
    for j=1:length(R)
        Lambda(i,j) = w(i)*2*pi*R(j)/(10*60); %Lambda at rated power
    end
end

lim=16/27;



%Plot results
figure()
for i=1:length(w)
    yyaxis left
    plot(R,Lambda(i,:))
    ylabel('Tip speed ratio \lambda [-]')
    xlabel('Blade radius [m]')
    legend('\Omega=450 rpm','\Omega=500 rpm','\Omega=550 rpm')

    yyaxis right
    plot(R,Cp)
    h = refline(0,lim);
    h.Color = 'black';
    h.LineStyle = '-.';
    hold on
    ylabel('Rotor power coefficient [-]')
    legend('\Omega=450 rpm','\Omega=500 rpm','\Omega=550 rpm')
    ylim([0 1])
end
grid on
