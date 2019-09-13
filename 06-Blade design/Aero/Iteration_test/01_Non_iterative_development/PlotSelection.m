function Plots = PlotSelection

%Function to select the plots to be made
%true to be printed
%false to be avoided

All_plots = 0;

%%Blade results as a function of lambda 
%Rotor orque and power vs omega, as a function of tip speed ratio
Plots.TorqueOmegaRotor = All_plots;

%%Operational points
%Curves resulting of the intersection between generator and rotor
%torque-omega curves
Plots.OperationalPoints = All_plots;

%%AEP
%Energy produced in one year probability curve
Plots.AEP = All_plots;

%%Chord and Twist
Plots.CT = All_plots;

end

