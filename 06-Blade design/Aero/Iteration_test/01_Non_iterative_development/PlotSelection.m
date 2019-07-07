function Plots = PlotSelection;

%Function to select the plots to be made
%true to be printed
%false to be avoided

%%Blade results as a function of lambda 
%Rotor orque and power vs omega, as a function of tip speed ratio
Plots.TorqueOmegaRotor = 0;

%%Operational points
%Curves resulting of the intersection between generator and rotor
%torque-omega curves
Plots.OperationalPoints = 0;

%%AEP
%Energy produced in one year probability curve
Plots.AEP = 0;

end

