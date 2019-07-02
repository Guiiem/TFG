function time_to_start = start_calc_opt(data,chord, twist, J, delt)

% Determine the two constants (c1, c2) in the equation for
% d(lambda)/dt

tw_rad = pi/180*twist; % Twist must be in radians for starting calcs

% constant terms in front of equation 4.8
c1 = 1.2*data.numb*data.u_start^2*data.r_tip^3*data.delr;
% constant terms in front of equation 4.9
c2 = data.r_tip/(J*data.u_start);
lambda = 0.0; time = 0.0;
Fn_3=0.0; Fn_2=0.0; Fn_1=0.0;

while (1) % Use Adams Moulton for integration
	fn_0 = deriv(chord, tw_rad, lambda, data);
	Fn_0 = c2*(c1*fn_0 - data.tor_resis);
    lam_pred=lambda + delt*(55*Fn_0-59*Fn_1+37*Fn_2-9*Fn_3)/24;
	fp_1 = deriv(chord, tw_rad, lam_pred, data);
	Fp_1 = c2*(c1*fp_1 - data.tor_resis);
    del_l = delt*(9*Fp_1+19*Fn_0-5*Fn_1+Fn_2)/24;
    if (lambda + del_l > data.lambda_start)
		time_to_start = time + delt*(data.lambda_start - lambda)/del_l;
		return
    end
    lambda = lambda + del_l;
	time = time + delt  ;
	if (time > data.t_max) % Penalise slow blades
		time_to_start = 5.5*data.t_max + 1.0 ;
		return
	end
	if (fp_1 < -1) 
		time_to_start = 4.5*data.t_max + 1.0;
		return
	end
	if (fp_1 < 0.0 | lambda < 0) 
		time_to_start = 3.5*data.t_max + 1.0;
		return
	end
	Fn_3=Fn_2;
    Fn_2=Fn_1;
    Fn_1=Fn_0;
end
end %start_calc  
