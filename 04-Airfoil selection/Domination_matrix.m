%Programa per fer la domination matrix i triar airfoil

Data = [0.67	0.99	0.44	0.33	0.73	0.76	0.48
0.65	0.55	0.66	0.68	0.93	0.48	0.49
0.75	0.84	0.34	0.61	0.90	0.48	0.46
0.88	1.00	0.00	0.04	0.65	0.48	0.53
0.47	0.78	0.29	0.86	1.00	0.86	0.41
0.49	0.69	0.69	0.91	1.00	0.71	0.38
0.43	0.92	0.39	0.90	0.95	1.00	0.37
0.95	0.58	0.71	0.00	0.72	0.57	0.76
0.85	0.51	0.64	0.23	0.99	0.58	0.83
0.56	0.36	0.86	1.00	0.66	0.33	0.37
0.78	0.60	0.96	0.53	0.92	0.43	0.51
0.74	0.54	0.98	0.52	0.91	0.49	0.52
0.65	0.49	1.00	0.72	0.93	0.42	0.46
0.70	0.57	0.51	0.92	0.87	0.24	0.59
0.86	0.62	0.32	0.25	0.36	0.43	0.48
0.68	0.73	0.43	0.75	0.56	0.52	0.39
1.00	0.64	0.25	0.06	0.00	0.27	0.52
0.66	0.49	0.78	0.80	0.73	0.42	1.00
];

Pes = [ 7 3 1 1 3 3 3 ]; %Pes de cada variable


dim = 18; %Numero d'alternatives que tinc

Resultat = zeros(dim,dim);

for i=1:dim
    for j=1:dim
        for k=1:7
            if(Data(i,k)>Data(j,k))
                Resultat(i,j)=Resultat(i,j)+(Data(i,k)-Data(j,k))*Pes(k);
            end
        end
    end    
end

Suma_fil = sum(Resultat.');
Suma_colm = sum(Resultat);

Valors = Suma_fil./Suma_colm;

[y,i] = max(Valors)




