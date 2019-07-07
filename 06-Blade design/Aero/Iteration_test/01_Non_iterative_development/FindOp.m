function Operation = FindOp(Results,Data,lambda,Plots)

Uin = Data.u_start;
Uout = Data.u_cut_out;

U = Uin:1:Uout;

Operation.Q = zeros(length(U),1);
Operation.Omega = zeros(length(U),1);
Operation.Lambda = zeros(length(U),1);


for j=1:length(U)
    
    %Torque-omega generator data
    Omega_G = Data.gen_omega;
    Q_G = Data.gen_torque;

    %Torque-omega rotor data
    Omega_R = Results.Omega(j,:);
    Q_R = Results.Q(j,:);
    
    %Initialize the variables that will test the results for each iteration
    OperationOmega_tst = 0;
    OperationQ_tst = 0;
    OperationLambda_tst = 0;
    res = 0; %Number of results (crossing points) for each iteration
    Slope_tst = 0; %Derivative of the solution in the crossing point

    for i=1:length(Omega_R)
        Q_G_new(i) = interp1(Omega_G,Q_G,Omega_R(i));
    end

    pos = 0;
    neg = 0;

    difpos = zeros(length(Omega_R),1);
    difneg = zeros(length(Omega_R),1);

    %Compute the torque difference at each omega
    for i=1:length(Omega_R)
        diff = Q_R(i)-Q_G_new(i);
        if diff>0
            pos = pos+1;
            difpos(i) = diff;
            difneg(i) = 1e10;
        end
        if diff<0
            neg = neg+1;
            difneg(i) = diff;
            difpos(i) = 1e10;
        end        
    end

    %If the torque is higher/lower at each omega, the curves do not interesect
    if(pos==0 || neg ==0)
      %  disp('No convergence at wind speed')
      %  disp(U(j))
        Operation.Q(j)=0;
        Operation.Omega(j)=0;
    end
    
    %New intersection method: evaluate at each section. If there is more
    %than one result, take the one with the rotor torque-omega with
    %negative slope
    
    if(pos~=0 && neg~=0)
        for i=1:(length(Omega_R)-1)
           
            Delta_Omega = Omega_R(i+1)-Omega_R(i);

            m1 = (Q_R(i+1)-Q_R(i))/Delta_Omega;
            m2 = (Q_G_new(i+1)-Q_G_new(i))/Delta_Omega;

            OperationOmega = (Q_G_new(i)-Q_R(i))/(m1-m2)+Omega_R(i);
            OperationQ = (OperationOmega-Omega_R(i))*m1+Q_R(i);
            
            %If it is in the range we are looking in, we have a solution 
            if (OperationOmega >= Omega_R(i) && OperationOmega <= Omega_R(i+1))
                
                res = res + 1;                
                OperationOmega_tst (res) = OperationOmega;
                OperationQ_tst(res) = OperationQ;
                Slope_tst(res) = m1;
                
                x = (Operation.Omega(j)-Omega_R(i))/Delta_Omega;
                OperationLambda_tst(res) = lambda(i)*x + lambda(i+1)*(1-x);
            end
        end
    end
    
    %Once we have obtained all the crossing points for two single
    %torque-omega curves, select the point with the negative slope (if more
    %than one result has been obtained)
    
    if (res == 1)
        
        Operation.Omega(j) = OperationOmega_tst;
        Operation.Q(j) = OperationQ_tst;
        Operation.Lambda(j) = OperationLambda_tst;
    end
    
    if (res > 1)
        for i=1:length(res)
            if (Slope_tst(i) <= 0)
                Operation.Omega(j) = OperationOmega_tst(res);
                Operation.Q(j) = OperationQ_tst(res);
                Operation.Lambda(j) = OperationLambda_tst(res);
            end
        end
    end
                

%     if(pos~=0 && neg~=0)
%         [M,I] = min(difpos(:));
%         [M,Y] = min(abs(difneg(:)));
%         i = min(I,Y);
%         Delta_Omega = Omega_R(i+1)-Omega_R(i);
%         
%         m1 = (Q_R(i+1)-Q_R(i))/Delta_Omega;
%         m2 = (Q_G_new(i+1)-Q_G_new(i))/Delta_Omega;
%         
%         Operation.Omega(j) = (Q_G_new(i)-Q_R(i))/(m1-m2)+Omega_R(i);
%         Operation.Q(j) = (Operation.Omega(j)-Omega_R(i))*m1+Q_R(i);
%         
%         %Find the operation lambda
%         x = (Operation.Omega(j)-Omega_R(i))/Delta_Omega;
%         Operation.Lambda(j) = lambda(i)*x + lambda(i+1)*(1-x);
%     end

end

%Finally compute power and power coefficent vs wind speed
Operation.P = Operation.Omega.*Operation.Q*2*pi/60;
Operation.Cp = Operation.P(:,1)'./(0.5*1.225*pi*Data.r_tip^2.*U(1,:).^3);

if Plots.OperationalPoints == true 
    figure('Name','Operation results with the generator curve')
    
    subplot(221)
    yyaxis left
    plot(Operation.Omega,Operation.Q)
    xlabel('Omega [rpm]')
    ylabel('Torque [Nm]')
    yyaxis right
    plot(Operation.Omega,Operation.Lambda)
    ylabel('Lambda [-]')

    subplot(222)
    yyaxis left
    plot(U,Operation.Q)
    hold on
    plot(U,Operation.Lambda)
    xlabel('Wind speed [m/s]')
    ylabel('Torque [Nm], Lambda[-]')
    yyaxis right
    plot(U,Operation.Omega)
    ylabel('Omega [rpm]')
    legend('Torque','Lambda')
       
    subplot(2,2,[3,4])
    yyaxis left
    plot(U,Operation.P)
    xlabel('Wind speed [m/s]')
    ylabel('Power [W]')
    yyaxis right
    plot(U,Operation.Cp);
    ylabel('Power coefficient [-]')
    
end
    
    

end
