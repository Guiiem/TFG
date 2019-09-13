function J = inertia(chord,twist,data)

%  Area  area of blade aerofoil (no units)
%  rho_blade density of blade material (kg m^-3)

Area = data.area; 
A2 = Area^2;
J1 = sum((chord.*data.rad).^2);
J2 = sum((cosd(twist).^2) .* (chord.^4)) + A2*sum((sind(twist).^2) .* (chord.^4));
% Convert to dimensioned torque and add drive train + generator contribution
J = data.numb*data.delr*data.rho_blade*Area*data.r_tip^5*(J1 + J2/12) + data.J_gen;
end