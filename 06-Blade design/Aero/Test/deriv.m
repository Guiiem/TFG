function dlamdt = deriv(chord, twist, lambda, data)
rad=data.rad;
lamr = lambda.*rad;
lamr2=lamr.*lamr;
tmp = sqrt(1 + lamr2).*sin(twist).*(cos(twist) - lamr.*sin(twist)).*...
		chord.*rad;
dlamdt=sum(tmp);