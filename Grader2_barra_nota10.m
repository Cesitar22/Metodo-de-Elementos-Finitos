clear all;
clc;
%% Datos de entrada
L = 1; a0 = L/40; E = 1000; P = 10; f0 =100;

%posibilidad de calcular los desplazamientos de manera simbólica
% syms x; %declara x como una variable simbólica
% f_x = f0*(x/L)^3;
% a_x = a0*exp(-x/(2*L));
% A_x = a_x*a_x;
% R = P + int(f_x, x, [0, L]);
% F_x = R - int(f_x, x, [0, x]); % ó hacer: sigma_x = P + int(f_x, x, [x, L]);
% sigma_x = F_x/A_x;
% eps_x = sigma_x/E;
% u_x = int(eps_x, x, [0, x]);
% u_fuerte = round(subs(u_x, x, L), 4);

%cálculo de la u_fuerte numéricamente
x_fuerte = L; %introducir el valor de x en el cuál se desea obtener la solución fuerte
R = P + f0*L/4;
A_x = a0^2*exp(-x_fuerte/L);
sigma_x = (R-(f0/(4*L^3)*x_fuerte^4))/(A_x);
eps_x = sigma_x/E;
u_x = 1/(E*a0^2)*(R*L*(exp(x_fuerte/L)-1)-f0/(4*L^3)*(24*L^5*exp(x_fuerte/L)-24*L^5-24*L^4*x_fuerte*exp(x_fuerte/L)+12*L^3*x_fuerte^2*exp(x_fuerte/L)-4*L^2*x_fuerte^3*exp(x_fuerte/L)+L*x_fuerte^4*exp(x_fuerte/L)));
u_fuerte = round(u_x, 4)

%% solucion MEF con elementos cuadraticos (3 barras de 3 nodos cada una); barra-viga
chi_pi1 = 0; w_pi1 = 2; %cuadratura n = 1
chi_pi2 = [-1 1]*1/sqrt(3); w_pi2 = [1 1]; %cuadratura n = 2
chi_pi3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_pi3 = [5/9 8/9 5/9]; %cuadratura n = 3

IPK = 3; %selección de la cantidad de puntos de integración por cada elemento, para la matriz de rigidez
if IPK == 1;
    chi_pik = chi_pi1; w_pik = w_pi1;
elseif IPK == 2;
    chi_pik = chi_pi2; w_pik = w_pi2;
else IPK == 3;
    chi_pik = chi_pi3; w_pik = w_pi3;
end

IPF = 3; %selección de la cantidad de puntos de integración por cada elemento, para el vector de fuerzas nodales equivalentes
if IPF == 1;
    chi_pif = chi_pi1; w_pif = w_pi1;
elseif IPF == 2;
    chi_pif = chi_pi2; w_pif = w_pi2;
else IPF == 3;
    chi_pif = chi_pi3; w_pif = w_pi3;
end

% Datos de malla
Le1 = L/3; Le2 = L/3; Le3 = L/3;
Je1 = Le1/2; Je2 = Le2/2; Je3 = Le2/2;
x_nod = [0 L/6 2*L/6 3*L/6 4*L/6 5*L/6 L]; %posición de los nodos en coordenadas globales
A_nodal1 = a0^2*[exp(-x_nod(1)/L); exp(-x_nod(2)/L); exp(-x_nod(3)/L)]; %área en los nodos del elemento 1
A_nodal2 = a0^2*[exp(-x_nod(3)/L); exp(-x_nod(4)/L); exp(-x_nod(5)/L)];
A_nodal3 = a0^2*[exp(-x_nod(5)/L); exp(-x_nod(6)/L); exp(-x_nod(7)/L)];
f_nodal1 = f0*[(x_nod(1)/L)^3; (x_nod(2)/L)^3; (x_nod(3)/L)^3]; %valor de la fuerza distribuida en los nodos del elemento 1
f_nodal2 = f0*[(x_nod(3)/L)^3; (x_nod(4)/L)^3; (x_nod(5)/L)^3];
f_nodal3 = f0*[(x_nod(5)/L)^3; (x_nod(6)/L)^3; (x_nod(7)/L)^3];

% Rigidez local de los elementos
% ELEMENTO 1
Ke1 = zeros(3);
for ip=1:IPK
    chi = chi_pik(ip); w = w_pik(ip);
    Be1 = (1/Je1)*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_chi = Ne*A_nodal1;
    Ke1 = Ke1 + Be1'*E*A_chi*Be1*Je1*w;
end
Ke1;

%ELEMENTO 2
Ke2 = zeros(3);
for ip=1:IPK
    chi = chi_pik(ip); w = w_pik(ip);
    Be2 = (1/Je2)*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_chi = Ne*A_nodal2;
    Ke2 = Ke2 + Be2'*E*A_chi*Be2*Je2*w;
end
Ke2;

%ELEMENTO 3
Ke3 = zeros(3);
for ip=1:IPK
    chi = chi_pik(ip); w = w_pik(ip);
    Be3 = (1/Je3)*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_chi = Ne*A_nodal3;
    Ke3 = Ke3 + Be3'*E*A_chi*Be3*Je3*w;
end
Ke3;

% Ensamblaje: Rigidez global
K_quest = zeros(7); %inicilización de la matriz 7x7 (pues hay 7 nodos en total) como una matriz de cerospues hay 7 nodos en total
gdl_e1 = [1 2 3];
K_quest(gdl_e1, gdl_e1) = Ke1 + K_quest(gdl_e1, gdl_e1);
gdl_e2 = [3 4 5];
K_quest(gdl_e2, gdl_e2) = Ke2 + K_quest(gdl_e2, gdl_e2);
gdl_e3 = [5 6 7];
K_quest(gdl_e3, gdl_e3) = Ke3 + K_quest(gdl_e3, gdl_e3);

K_quest;

% fuerza nodal local de los elementos
% ELEMENTO 1
Fe1 = zeros(3, 1); %inicialización del vector de fuerzas nodales equivalentes  del elemento 1
for ip=1:IPF
    chi = chi_pif(ip); w = w_pif(ip);
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    f_chi = Ne*f_nodal1; %fuerza externa virgulilla (en formulación isoparaméttrica)
    Fe1 = Fe1 + Ne'*f_chi*Je1*w;
end
Fe1; %fuerzas externas nodales equivalentes del elemento 1

% ELEMENTO 2
Fe2 = zeros(3, 1); %inicialización del vector de fuerzas nodales equivalentes  del elemento 2
for ip=1:IPF
    chi = chi_pif(ip); w = w_pif(ip);
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    f_chi = Ne*f_nodal2; %fuerza externa virgulilla (en formulación isoparaméttrica)
    Fe2 = Fe2 + Ne'*f_chi*Je2*w;
end
Fe2; %fuerzas externas nodales equivalentes del elemento 2

% ELEMENTO 3
Fe3 = zeros(3, 1); %inicialización del vector de fuerzas nodales equivalentes  del elemento 3
for ip=1:IPF
    chi = chi_pif(ip); w = w_pif(ip);
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    f_chi = Ne*f_nodal3; %fuerza externa virgulilla (en formulación isoparaméttrica)
    Fe3 = Fe3 + Ne'*f_chi*Je3*w;
end
Fe3; %fuerzas externas nodales equivalentes del elemento 3

% fuerza nodal global (emsamblaje)
f_quest = zeros(7, 1);
f_quest(gdl_e1, 1) = Fe1 + f_quest(gdl_e1, 1);
f_quest(gdl_e2, 1) = Fe2 + f_quest(gdl_e2, 1);
f_quest(gdl_e3, 1) = Fe3 + f_quest(gdl_e3, 1);

f_direct = zeros(7, 1); %fuerzas directamente aplicadas en los nodos
f_direct(7) = P;

f_quest = f_quest + f_direct;

% Condiciones de contorno: calcula desplazamientos nodales
gdl_T = [1:7]; %gdl totales
gdl_R = 1; %gdl restringidos
gdl_L = setdiff(gdl_T, gdl_R); %gdl libres
KLL = K_quest(gdl_L, gdl_L); %matriz de rigidez reducida
f_L = f_quest(gdl_L);
uL_quest = KLL\f_L; %vector de desplazamientos en nodos libres
U = zeros(length(gdl_T),1);%vector de desplazamientos completo
U(gdl_L) = uL_quest;

%display de RESULTADOS
K_quest = round(K_quest, 2)
f_quest = round(f_quest, 2)
uL_quest = round(uL_quest, 2)

% El error relativo % en desplazamientos en el extremo libre entre la solución fuerte y la solución mef (error_u)
error_u = round(100*(abs(u_fuerte-U(7)))/u_fuerte, 2) %en tanto por 100


%% Determinar la deformacion y tension

%% Tensión máxima de la solución fuerte
x_n = [0:0.001:L];
sigma_xi = zeros(length(x_n), 1);
for n = 1:length(x_n)
    x_i = x_n(n);
    A_xi = a0^2*exp(-x_i/L);
    sigma_xi(n) = (R-(f0/(4*L^3)*x_i^4))/(A_xi);
end
sigma_xi;
sigmamax_fuerte = max(sigma_xi)

%% Deformación y tensión de la solución débil
coord_nodal = [0 L/6 2*L/6 3*L/6 4*L/6 5*L/6 L]; % coordenadas de los nodos
%ELEMENTO 1
gdle1 = [1; 2; 3];
x_e1 = coord_nodal(gdle1); % coordenadas nodales del elemento 1
u_e1 = U(gdle1); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en elemento 1

% deformaciones y tensiones en los puntos de integración
epsilon_ip_e1 = zeros(2,1); %deformación del elemento 1 en los puntos de integración de K
for ip = 1:2
    chi = chi_pi2(ip);%datos de integración    
    Be1 = Je1^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_ip_e1(ip) = Be1*u_e1;    
end
% deformaciones y tensiones en los nodos
epsilon_nod_e1 = zeros(3,1); %deformación del elemento 1 en los puntos de integración de K
chi_nodal = [-1; 0; 1];
for n = 1:3
    chi = chi_nodal(n);%datos de integración    
    Be1 = Je1^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_nod_e1(n) = Be1*u_e1;
end
simga_nod_e1 = epsilon_nod_e1*E; %tensión del elemento 1 en los puntos de integración de K

%ELEMENTO 2
gdle2 = [3; 4; 5];
x_e2 = coord_nodal(gdle2); % coordenadas nodales del elemento 2
u_e2 = U(gdle2); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en elemento 2

% deformaciones y tensiones en los puntos de integración
epsilon_ip_e2 = zeros(2,1); %deformación del elemento 2 en los puntos de integración de K
for ip = 1:2
    chi = chi_pi2(ip);%datos de integración    
    Be2 = Je2^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_ip_e2(ip) = Be2*u_e2;    
end
% deformaciones y tensiones en los nodos
epsilon_nod_2 = zeros(3,1); %deformación del elemento 2 en los puntos de integración de K
chi_nodal = [-1; 0; 1];
for n = 1:3
    chi = chi_nodal(n);%datos de integración    
    Be2 = Je2^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_nod_e2(n) = Be2*u_e2;
end
simga_nod_e2 = (epsilon_nod_e2*E)'; %tensión del elemento 2 en los puntos de integración de K


%ELEMENTO 3
gdle3 = [5; 6; 7];
x_e3 = coord_nodal(gdle3); % coordenadas nodales del elemento 1
u_e3 = U(gdle3); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en elemento 1

% deformaciones y tensiones en los puntos de integración
epsilon_ip_e3 = zeros(2,1); %deformación del elemento 1 en los puntos de integración de K
for ip = 1:2
    chi = chi_pi2(ip);%datos de integración    
    Be3 = Je3^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_ip_e3(ip) = Be3*u_e3;    
end
% deformaciones y tensiones en los nodos
epsilon_nod_e3 = zeros(3,1); %deformación del elemento 1 en los puntos de integración de K
chi_nodal = [-1; 0; 1];
for n = 1:3
    chi = chi_nodal(n);%datos de integración    
    Be3 = Je3^-1*[chi-1/2 -2*chi chi+1/2]; %matriz cinemática en formato numerico
    epsilon_nod_e3(n) = Be3*u_e3;
end
simga_nod_e3 = epsilon_nod_e3*E; %tensión del elemento 1 en los puntos de integración de K

sigmamax_nod = max(max([simga_nod_e1 simga_nod_e2 simga_nod_e3])) %máximo de sigma_ip_e1

err_sigmamax_nod = round(100*(abs(sigmamax_nod-sigmamax_fuerte))/sigmamax_fuerte, 2)

%% Representa la tensión nodal a lo largo de la barra
figure(1)
grafstressnod = plot(x_e1,simga_nod_e1,'go-','LineWidth',2) %representa el elemento 1 y almacena en la variable grafstressnod
hold on %mantiene los datos de la gráfica anterior cuando representes los siguientes elementos
grafstressnod = plot(x_e2,simga_nod_e2,'go-','LineWidth',2) %representa el elemento 2
grafstressnod = plot(x_e3,simga_nod_e3,'go-','LineWidth',2) %representa el elemento 3
title('Tension en solución fuerte vs MEF en nodos','Fontsize',10);
% Representa la tensión de solución fuerte a lo largo de la barra
x_ = 0:0.1:1; %puntos soporte para la representación fuerte
for i=1:length(x_)
    x = x_(i);
    ufuerte(i) = 1/(E*a0^2)*(R*L*(exp(x/L)-1)-f0/(4*L^3)*(24*L^5*exp(x/L)-24*L^5-24*L^4*x*exp(x/L)+12*L^3*x^2*exp(x/L)-4*L^2*x^3*exp(x/L)+L*x^4*exp(x/L)));
    epsfuerte(i) = (R-(f0/(4*L^3)*x^4))/(E*(a0^2*exp(-x/L)));
    sigmafuerte(i) =  E*epsfuerte(i); 
end
grafstress_strong = plot(x_,sigmafuerte,'r-','LineWidth',2)

err_sigmamax_nod = err_sigmamax_nod; %error relativo entre el máximo nodal y el máximo de la solución fuerte

% Representa el desplazamiento a lo largo de la barra
figure(2)
graf_u_strong = plot(x_,ufuerte,'r-','LineWidth',2) %representa el desplazamiento fuerte frente a la longitud L
hold on
graf_u_mef = plot(coord_nodal,U,'go','LineWidth',2) %representa las coordenadas de todos los nodos frente al desplazamiento mef U
title('Desplazamiento fuerte vs MEF','Fontsize',10);
