function AEP = AnnualEnergyProduction(Operation,Data,Plots)

%Generation of the wind speed vector
U = Operation.U;

%Weibull parameters
A = Data.Weibull_A;
k = Data.Weibull_k;

%Compute the probability of the Weibull distribution
Prob = k/A.*(U/A).^(k-1).*exp(-(U/A).^k);

%Express the probability in terms of hours per year
Hours = Prob*24*365;

%Energy produced for each wind speed in one year [W*h]
E_prod = Hours.'.*Operation.P;

%Numerical integration of the energy production probability
AEP = trapz(U,E_prod);

if Plots.AEP == true
    
    figure('Name','Annual Energy production')
    yyaxis left
    plot(U,E_prod/1000)
    xlabel('Wind speed [m/s]')
    ylabel('Energy produced in one year [kW*h]')
    yyaxis right
    plot(U,Operation.P)
    ylabel('Power output [W]')
    
end







end