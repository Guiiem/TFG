function Results = Curves(Cp,Cq,Data,lambda,Plots);

%First, plot the power and torque curves vs wind speed as function of
%lambda

Uin = Data.u_start;
Uout = Data.u_cut_out;
Rho = 1.225;
R = Data.r_tip;

U = Uin:1:Uout;

for i=1:length(U)
    for j=1:length(lambda)
        Results.P(i,j) = Cp(j)*0.5*Rho*pi*(R)^2*U(i)^3; %Aerodynamic power for each (wind speed,lambda)
        Results.Omega(i,j) = lambda(j)*U(i)/(R); %Rotor speed for each (wind speed, lambda)
        Results.Q(i,j) = Results.P(i,j)/Results.Omega(i,j); %Aerodynamic torque for each (wind speed, lambda)
    end
end


if(Plots.TorqueOmegaRotor == true)

    figure('Name','Blade results as a function of lambda')
    subplot(221)
    plot(lambda,Cp)
    xlabel('Lambda')
    ylabel('Power coefficent')

    subplot(222)
    plot(lambda,Cq)
    xlabel('Lambda')
    ylabel('Torque coefficent')

    subplot(2,2,3)
    title('Torque vs Omega for each wind speed');
    for i=1:(length(U)/2)
        plot(Results.Omega(i*2-1,:),Results.Q(i*2-1,:));
        hold on
        legendInfo{i}=['U: ' num2str(U(i*2-1))];
    end
    hold on
    plot(Data.gen_omega,Data.gen_torque)
    legend(legendInfo,'Generator curve')
    xlabel('Omega [rpm]')
    ylabel('Torque [Nm]')

    subplot(2,2,4)
    for i=1:(length(U)/2)
        plot(Results.Omega(i*2-1,:),Results.P(i*2-1,:));
        hold on
        legendInfo{i}=['U: ' num2str(U(i*2-1))];
    end
    hold on
    plot(Data.gen_omega,Data.gen_torque.*Data.gen_omega*2*pi/(60))
    legend(legendInfo,'Generator curve')
    xlabel('Omega [rpm]')
    ylabel('Power [W]')
    
end

end

