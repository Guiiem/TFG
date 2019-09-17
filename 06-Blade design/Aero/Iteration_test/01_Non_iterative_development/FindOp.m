function Operation = FindOp(Results,Data,lambda,Plots,Power_Curve)

Uin = Data.u_start;
Uout = Data.u_cut_out;
Ustep = Data.u_step;

U = Uin:Ustep:Uout;

Operation.Q = zeros(length(U),1);
Operation.Omega = zeros(length(U),1);
Operation.Lambda = zeros(length(U),1);
Operation.T = zeros(length(U),1);



for j=1:length(U)
    
    %Torque-omega generator data
    Omega_G = Data.gen_omega;
    Q_G = Data.gen_torque;

    %Torque-omega rotor data
    Omega_R = Results.Omega(j,:);
    Q_R = Results.Q(j,:);
    
    %Thrust curve
    Thrust = Results.T(j,:);
    
    %Initialize the variables that will test the results for each iteration
    OperationOmega_tst = 0;
    OperationQ_tst = 0;
    OperationLambda_tst = 0;
    Thrust_tst = 0;
    res = 0; %Number of results (crossing points) for each iteration
    Slope_tst = 0; %Derivative of the solution in the crossing point

    for i=1:length(Omega_R)
        %Generator and rotor omega curves should have the same length
        i_up = find (Omega_G >= Omega_R(i));
        
        if (~isempty(i_up))            
            factor = (Omega_R(i) - Omega_G(i_up-1))/(Omega_G(i_up)-Omega_G(i_up-1));
            Q_G_new(i) = Q_G(i_up)*factor + Q_G(i_up-1)*(1-factor);
            i_last = i_up;
        end
        
        %if the rotor curve is outside the region
        if (isempty(i_up))
%             m = (Q_G(i_last-1)-Q_G(i_last))/(Omega_G(i_last)-Omega_G(i_last-1));
%             Q_G_new(i) = Q_G(i_last) + (Omega_R(i)-Omega_G(i_last))*m;
              Q_G_new(i) = NaN;
        end
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
                
                x = (OperationOmega-Omega_R(i))/Delta_Omega;
                OperationLambda_tst(res) = lambda(i)*(1-x) + lambda(i+1)*(x);
                
                %Also compute the Thrust in this region
                Thrust_tst(res) = Thrust(i)*(1-x) + Thrust(i+1)*x;
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
        Operation.T(j) = Thrust_tst;    
        
    elseif (res > 1)
        for i=1:length(res)
            if (Slope_tst(i) <= 0)
                Operation.Omega(j) = OperationOmega_tst(res);
                Operation.Q(j) = OperationQ_tst(res);
                Operation.Lambda(j) = OperationLambda_tst(res);
                Operation.T(j) = Thrust_tst(res);
                
            end
        end
        
        %If there is no intersection and the wind speed is high, set the
        %current wind speed as the new cut out wind speed.
        
    elseif (res == 0 && U(j)>=8)
            U = Uin:Ustep:U(j-1);
            Operation.U = U;
            Operation.Omega = Operation.Omega(1:j-1);
            Operation.Lambda = Operation.Lambda(1:j-1);
            Operation.Q = Operation.Q(1:j-1);
            Operation.T = Operation.T(1:j-1);
            
            %U(j) = Uout; %End loop
            break
            
    elseif (res == 0 && U(j)<8)
        Operation.Omega(j) = NaN;
        Operation.Q(j) = NaN;
        Operation.Lambda(j) = NaN;
        Operation.T(j) = NaN;          
        
    else
        disp('Passa algo');
    end

end

%Then compute power and power coefficent vs wind speed
Operation.P = Operation.Omega.*Operation.Q*2*pi/60;
Operation.Cp = Operation.P(:,1)'./(0.5*1.225*pi*Data.r_tip^2.*U(1,:).^3);

%Compute thrust coefficient
Operation.Ct = Operation.T(:,1)'./(0.5*1.225*pi*Data.r_tip^2.*U(1,:).^2);



%We know the lambda at each operational point. Compute the mean angle of
%attack, lift and drag coefficient there
for i=1:length(Operation.Lambda)
    lam = Operation.Lambda(i);
    if ~isnan(lam)
        Del_lam = Power_Curve.Lambda(2)-Power_Curve.Lambda(1);
        %Find the index of the power curve lambda vector
        [~,ind] = find(Power_Curve.Lambda > lam, 1);
        if ind == 1
            Operation.MeanAoA(i) = Power_Curve.MeanAoA(ind);
            Operation.MeanCl(i) = Power_Curve.MeanCl(ind);
        else
            x = (lam-Power_Curve.Lambda(ind-1))/Del_lam;
            Operation.MeanAoA(i) = Power_Curve.MeanAoA(ind-1)*(1-x) + Power_Curve.MeanAoA(ind)*(x);
            Operation.MeanCl(i) = Power_Curve.MeanCl(ind-1)*(1-x) + Power_Curve.MeanCl(ind)*(x);
        end
    else
        Operation.MeanAoA(i) = NaN;
        Operation.MeanCl(i) = NaN;
    end
end

%Eliminate the values of U before cut-in
for i = 1:length(U)
    if isnan(Operation.Lambda(i))
       ind = i;
    end 
end

fin = length(U);
U = U(ind+1:fin);
Operation.U = U;
Operation.Omega = Operation.Omega(ind+1:fin);
Operation.Lambda = Operation.Lambda(ind+1:fin);
Operation.P = Operation.P(ind+1:fin);
Operation.Q = Operation.Q(ind+1:fin);
Operation.T = Operation.T(ind+1:fin);
Operation.Cp = Operation.Cp(ind+1:fin);
Operation.Ct = Operation.Ct(ind+1:fin);
Operation.MeanAoA = Operation.MeanAoA(ind+1:fin);
Operation.MeanCl = Operation.MeanCl(ind+1:fin);









%Plots
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
       
    subplot(2,2,3)
    yyaxis left
    plot(U,Operation.P)
    xlabel('Wind speed [m/s]')
    ylabel('Power [W]')
    yyaxis right
    plot(U,Operation.Cp);
    ylabel('Power coefficient [-]')
    
    subplot(2,2,4)
    yyaxis left
    plot(U,Operation.T)
    xlabel('Wind speed [m/s]')
    ylabel('Trhust [N]')
    yyaxis right
    plot(U,Operation.Ct);
    ylabel('Thrust coefficient [-]')
    
end
    
    

end
