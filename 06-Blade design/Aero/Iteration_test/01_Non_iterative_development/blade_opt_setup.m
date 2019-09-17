function [data]=blade_opt_setup;
% This function sets up the data structure "data"

% PARAMETERs for the evolutionary optimisation
% data.no_indiv  number of members of each generation, should be 1000 or more
data.no_indiv=10^2;
% data.no_gen    number of generation to evolve over
data.no_gen=1;
% data.crossover crossover probability for breeding new blades
data.crossover=0.1;
% data.max_age   age at which blades "die", i.e. are removed from the population
data.max_age=20;
% Aerofoil parameters
%  data.re_array  array holding the Reynolds numbers.
data.re_array = [120000 160000 200000];
%  data.num_re    number of values of Reynolds number in data file.
data.num_re=length(data.re_array);
%  data.re_log    array holding the log of the Reynolds numbers
data.re_log=log(data.re_array);
%  data.alpha_array   array holding the angles of attack.
aerofoil_in = load('9_S1223.in'); 
%  data.num_alpha number of values of angle of attack at each Re.
data.num_alpha = length(aerofoil_in); 
%  data.cl_array, data.cd_array    arrays holding lift and drag coefficients.
data.alpha =zeros(data.num_alpha,1);
data.cl_array=zeros(data.num_alpha,data.num_re);
data.cd_array=zeros(data.num_alpha,data.num_re);
for i=1:data.num_alpha %For each angle of attack
        data.alpha(i)=aerofoil_in(i,1);
        for j=1:data.num_re
            data.cl_array(i,j)=aerofoil_in(i,j*2);
            data.cd_array(i,j)=aerofoil_in(i,j*2+1);
        end
end
% data.del_alpha angle between successive Cl and Cd values
data.del_alpha=data.alpha(2)-data.alpha(1);
% data.numb      number of blades
data.numb=3;
% data.area     (aerofoil cross-sectional area)/c^2
data.area=659.89/(100^2);
% data.rho_blade     blade density (kg/m^3)
data.rho_blade=430;
% data.no_genes      number of genes (blade elements) per blade
data.no_genes=15;
% data.u_start   wind speed (m/s) for the starting calculation
data.u_start=2.0;
% data.u_step   wind speed steps (m/s) for all the calculations
data.u_step = 0.2;
% data.u_cut_out    wind speed (m/s) for cut out
data.u_cut_out = 16;
% data.u_power   wind sTpeed (m/s) for the power extraction calculation
data.u_power=10;
% data.lambda_start  the tip speed ratio (TSR) to "complete" starting
data.lambda_start=1;
% data.Cp_min allowable Cp
data.Cp_min=0.01;
% data.t_max maximum time allowed for starting, determined in program
data.t_max=600;
% data.lambda_power  TSR for power extraction (will be defined again later)
data.lambda_power=4.608;

% data.lambda_range     range of lambda to be studied
data.lambda_range = 0.7:0.02:4;

% data.r_tip     blade tip radius (m)
data.r_tip=0.5;
% data.r_hub     Hub radius
data.r_hub=0.1;
% data.max_chord     maximum chord on blade
data.max_chord=0.3;
% data.min_chord     minimum chord on blade
data.min_chord=0.03;
% data.max_twist     maximum blade twist (in degrees)
data.max_twist=50;
%data. min_twist     minimum blade twist (in degrees)
data.min_twist=-20;
% data.J_gen     contribution to moment of inertia from generator & drive train
data.J_gen=0.006;
% data.tor_resis resistive torque of generator and drive train (Nm)
data.tor_resis=0.0;

%Gearbox ratio
Gr = 1;

% datal.delt    timestep for strating calculation
data.delt = 0.15;

%data.gen_power    generator power curve
data.gen_power = [0 500];
% data.gen_torque   generator torque curve
data.gen_torque =  [0.42 500/(600*2*pi/60)*Gr];
% data.gen_omega    generator rpm curve
data.gen_omega = [0 600/Gr];
% data.Weibull_A     Weibull parameter A
data.Weibull_A = 4;
% data.Weibull_k     Weibull paralameter k
data.Weibull_k = 2;

% data.a_power_in    factor (<= 1) to scale fitness score for power
data.a_power_in=1;
% data.non_dom   index of non dominated blades
data.non_dom=1;
% data.unfit     index of unfit blades
data.unfit=1;
% data.dominated     index of dominated blades
data.dominated=2;
% data.no_nondom     number of non dominated blades
data.no_nondom=0;
% data.no_unfit      number of unfit blades
data.no_unfit=0;
% data.no_dom    number of dominated blades
data.no_dom=0;
% data.delr 	difference between consecutive blade elements
data.delr=(data.r_tip - data.r_hub)/data.no_genes;
% data.rad	radius of each blade element
data.rad=data.r_hub+0.5*data.delr:data.delr:data.r_tip-0.5*data.delr;
data.rad=data.rad./data.r_tip;
%data.rad_boundary     radius of the boundary of each blade element [m]
data.rad_boundary = data.r_hub:data.delr:data.r_tip;
data.maxchord = zeros(1,data.no_genes); data.maxtwist = zeros(1,data.no_genes)
% data.fitfile 		String containing the file name to output the fittest blades to.
data.fitfile='fit_blade_out.dat';
%