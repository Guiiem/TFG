function Results = Curves(Cp,Cq,Ct,Data,lambda,Plots)

%Plot the power and torque curves vs wind speed as function of lambda

Uin = Data.u_start;
Uout = Data.u_cut_out;
Ustep = Data.u_step;
Rho = 1.225;
R = Data.r_tip;

U = Uin:Ustep:Uout;

for i=1:length(U)
    for j=1:length(lambda)
        Results.P(i,j) = Cp(j)*0.5*Rho*pi*(R)^2*U(i)^3; %Aerodynamic power for each (wind speed,lambda)
        Results.T(i,j) = Ct(j)*0.5*Rho*pi*R^2*U(i)^2; %Thrust for each (wind speed, lambda)
        Results.Omega(i,j) = (lambda(j)*U(i)/R)*60/(2*pi); %Rotor speed for each (wind speed, lambda) [rpm]
        Results.Q(i,j) = Results.P(i,j)/(Results.Omega(i,j)*2*pi/60); %Aerodynamic torque for each (wind speed, lambda)
    end
end


if(Plots.TorqueOmegaRotor == true)

    figure('Name','Blade results as a function of lambda')
    subplot(221)
    plot(lambda,Cp)
        xlim([0 6])
    xlabel('Lambda [-]')
    ylabel('Power coefficent [-]')

    subplot(222)
    yyaxis right
    plot(lambda,Cq)
            xlim([0 6])
    xlabel('Lambda [-]')
    ylabel('Torque coefficent [-]')
    yyaxis left
            xlim([0 6])
    plot(lambda,Ct)
    ylabel('Thrust coefficient [-]')

    subplot(2,2,4)
    title('Torque vs Omega for each wind speed');
    for i=1:(length(U)/10)
        plot(Results.Omega(i*10-9,:),Results.Q(i*10-9,:));
    %for i=1:(length(U))
     %   plot(Results.Omega(i,:),Results.Q(i,:));
        hold on
%        legendInfo{i}=['U: ' num2str(U(i*2-1))];
    end
    hold on
    xlim([0 800])
    ylim([0 15])
    plot(Data.gen_omega,Data.gen_torque)
   % legend(legendInfo,'Generator curve')
    xlabel('Omega [rpm]')
    ylabel('Torque [Nm]')

    Gen_Omega = Data.gen_omega(1):1:Data.gen_omega(2);
    Gen_Power = Gen_Omega.^2.*Data.gen_torque(2)./Data.gen_omega(2)*2*pi/60;
    
    
    subplot(2,2,3)
   for i=1:(length(U)/10)
        plot(Results.Omega(i*10-9,:),Results.P(i*10-9,:));
  %  for i=1:(length(U))
   %     plot(Results.Omega(i,:),Results.P(i,:));
   
        hold on
        legendInfo{i}=['U = ' num2str(U(i*10-9)) ' m/s'];
    end 
    hold on
    xlim([0 800])
    ylim([0 650])
    plot(Gen_Omega,Gen_Power)
    legend(legendInfo,'Generator')
    xlabel('Omega [rpm]')
    ylabel('Power [W]')
    x0=10;
    y0=10;
    width=700;
    height=400;
    set(gcf,'position',[x0,y0,width,height])

    
    
    figure()
    title('Torque vs Omega for each wind speed');
    for i=1:(length(U)/3)
        plot(Results.Omega(3*i-2,:),Results.Q(3*i-2,:));
        hold on
       legendInfo{i}=['U = ' num2str(U(i*3-2)) ' m/s'];
    end
    hold on
    plot(Data.gen_omega,Data.gen_torque)
    legend(legendInfo,'Generator curve')
    xlabel('Omega [rpm]')
    ylabel('Torque [Nm]')
    xlim([0 600])
    ylim([0 10])
    
end

end

