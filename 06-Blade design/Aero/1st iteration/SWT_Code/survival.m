function fitness= survival(data,bladechord,bladetwist)
delt = min(data.t_max/50, 0.25);
for i=1:data.no_indiv;
    cp(i)=power_calc_opt(data,bladechord(i,:),bladetwist(i,:),0);
    J = inertia(bladechord(i,:),bladetwist(i,:),data);
    Ts(i)=start_calc_opt(data,bladechord(i,:),bladetwist(i,:),J, delt);
end
fitness=[cp' Ts'];