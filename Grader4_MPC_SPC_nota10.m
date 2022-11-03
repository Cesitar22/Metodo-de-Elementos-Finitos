close all
clear all
clc

%% Datos de entrada
ne1 = 5; %número de elementos de la viga principal (Timoshenko)
ne2 = 2; %número de elementos de la riostra
gdln = 3;%grados de libertad por nodo axil, cortante y flector
nod_e = 2; % número de nodos de cada elemento (viga de dos nodos)

E = 70e9; nu = 0.3; L = 5; G = E/(2*(1+nu));%datos comunes a la viga principal y a la riostra
a1 = L/40; b1 = a1/2; t1 = a1/10; coef1 = 2.5; Lv = L; Le1 = Lv/ne1; Je1 = Le1/2;%viga principal (Timoshenko)
a2 = a1/5; t2 = a2/10; coef2 = 2; Lr = sqrt((L/2)^2+(L/3)^2); Le2 = Lr/ne2; Je2 = Le2/2;%riostra
alfar = atan((L/2)/(L/3)); %ángulo de inclinación de la riostra respecto a coordenadas globales
u5 = -L/100; %desplazamiento vertical hacia abajo impuesto en el nodo 5
x_MPC = L/3; %punto de aplicación de la MPC
k = 10^6; %parámetro de penalización

IPK = 3; %selección de la cantidad de puntos de integración por cada elemento, para la matriz de rigidez
IPF = 3; %selección de la cantidad de puntos de integración por cada elemento, para el vector de fuerzas nodales equivalentes
IPKr = 1; %IPK para integración reducida de la matriz de flexión y cortadura


%% Definición de variables
nn1 = ne1 + 1; %número de nodos de la viga principal
nn2 = ne2 + 1; %número de nodos de la riostra
nnt = nn1 + nn2; %número de nodos totales
nod_1 = [1:1:nn1]; %nodos que le corresponden a la viga principal
nod_2 =[nn1+1:1:nn1+nn2]; %nodos que le corresponden a la riostra
gdl_T = [1:nnt*gdln]; %gdl totales
% gdl_R = [1, 2, 3]; %gdl restringidos
% gdl_L = setdiff(gdl_T, gdl_R); %gdl libres
Rr = [cos(alfar) sin(alfar) 0; -sin(alfar) cos(alfar) 0; 0 0 1]; %matriz de rotación del elemento viga
Rr_ = [Rr zeros(3); zeros(3) Rr]; %matriz de rotación del elemento viga ampliada

%% Cálculos preliminares

Av = 2*t1*b1 + (a1-2*t1)*t1; %área de las secciones de la viga principal
Ar = a2^2 - (a2-2*t2)^2; %área de las secciones de la riostra
Av_c = Av/coef1; %área equivalente en energía para la viga principal (área efectiva en cortante)
Ar_c = Ar/coef2; %área equivalente en energía para la riostra (área efectiva en cortante)
Iv = (1/12)*t1*(a1-2*t1)^3 + 2*((1/12)*t1^3*b1 + t1*b1*(a1/2-t1/2)^2); %momento de inercia de la viga principal
Ir = (1/12)*(a2^4 - (a2 - 2*t2)^4); %momento de inercia de la riostra

A_nodalv = zeros(nn1, 1);
for i = 1:nn1
    A_nodalv(i,1) = Av;
end
A_nodalv; %valor del área en los nodos de la viga principal

I_nodalv = zeros(nn1, 1);
for i = 1:nn1
    I_nodalv(i,1) = Iv;
end
I_nodalv; %valor del área en los nodos de la viga principal

A_nodalr = zeros(nn2, 1);
for i = 1:nn2
    A_nodalr(i,1) = Ar;
end
A_nodalr; %valor del área en los nodos de la riostra

I_nodalr = zeros(nn2, 1);
for i = 1:nn2
    I_nodalr(i,1) = Ir;
end
I_nodalr; %valor del área en los nodos de la riostra

Av_c_nodal = zeros(nn1, 1);
for i = 1:nn1
    Av_c_nodal(i,1) = Av_c;
end
Av_c_nodal; %valor del área CORREGIDA en los nodos de la viga principal

Ar_c_nodal = zeros(nn2, 1);
for i = 1:nn2
    Ar_c_nodal(i,1) = Ar_c;
end
Ar_c_nodal; %valor del área CORREGIDA en los nodos de la riostra

gdl_nodalv = zeros(ne1,gdln*2); 
for i = 1:ne1-1
   for j = 1:gdln*2
           gdl_nodalv(1,j) = j; 
           gdl_nodalv(i+1,j) = gdl_nodalv(1,j)+3*i;
    end
end
gdl_nodalv; %grados de libertad que le corresponden a cada a ELEMENTO de la viga principal

gdl_nodalr = zeros(ne2,gdln*2); 
for i = 1:ne2-1
   for j = 1:gdln*2
           gdl_nodalr(1,j) = j + gdl_nodalv(ne1,gdln*2); 
           gdl_nodalr(i+1,j) = gdl_nodalr(1,j) + 3*i;
    end
end
gdl_nodalr; %grados de libertad que le corresponden a cada a ELEMENTO de la riostra

nodos_ev = zeros(ne1, nod_e);
for i = 1:ne1-1
    for j = 1:nod_e
        nodos_ev(1,j) = j;
        nodos_ev(i+1,j) = nodos_ev(1,j)+i;
    end
end
nodos_ev; %nodos de cada elemento de la viga principal(una fila es un elemento)

nodos_er = zeros(ne2, nod_e);
for i = 1:ne2-1
    for j = 1:nod_e
        nodos_er(1,j) = j;
        nodos_er(i+1,j) = nodos_er(1,j)+i;
    end
end
nodos_er; %nodos de cada elemento de la riostra(una fila es un elemento)


%% Cuadraturas de integración
chi_pi1 = 0; w_pi1 = 2; %cuadratura n = 1
chi_pi2 = [-1 1]*1/sqrt(3); w_pi2 = [1 1]; %cuadratura n = 2
chi_pi3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_pi3 = [5/9 8/9 5/9]; %cuadratura n = 3
chi_pi4 = [-sqrt((3-2*sqrt(6/5))/7) -sqrt((3+2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7)]; w_pi4 = [(18+sqrt(30))/36 (18-sqrt(30))/36 (18-sqrt(30))/36 (18+sqrt(30))/36]; %cuadratura n = 4

if IPK == 1;
    chi_pik = chi_pi1; w_pik = w_pi1;
elseif IPK == 2;
    chi_pik = chi_pi2; w_pik = w_pi2;
elseif IPK == 3;
    chi_pik = chi_pi3; w_pik = w_pi3;
else IPK == 4;
    chi_pik = chi_pi4; w_pik = w_pi4;
end

if IPKr == 1;
    chi_pikr = chi_pi1; w_pikr = w_pi1;
elseif IPKr == 2;
    chi_pikr = chi_pi2; w_pikr = w_pi2;
elseif IPKr == 3;
    chi_pikr = chi_pi3; w_pikr = w_pi3;
else IPKr == 4;
    chi_pikr = chi_pi4; w_pikr = w_pi4;
end

if IPF == 1;
    chi_pif = chi_pi1; w_pif = w_pi1;
elseif IPF == 2;
    chi_pif = chi_pi2; w_pif = w_pi2;
elseif IPF == 3;
    chi_pif = chi_pi3; w_pif = w_pi3;
else IPF == 4;
    chi_pif = chi_pi4; w_pif = w_pi4;
end

%% Cálculo de matriz de rigidez de la estructura K_global
%Emplea la matriz de rigidez del elemento viga de Timoshenko que aparece en las diapositivas 62 o 65

K_quest = zeros(length(gdl_T));
F = zeros(length(gdl_T),1);
gdl_flexion = [2, 3, 5, 6];
gdl_axil = [1 4];

%viga principal
for g = 1:ne1
    Ke = zeros(nod_e*gdln); %inicialización de la matriz de rigidez en axil y flexión de cada elemento
    Kef = zeros(nod_e*2); %inicialización de la matriz de rigidez en flexión de cada elemento
    Kec = zeros(nod_e*2);
    Kefc = zeros(nod_e*2);
    Kea = zeros(nod_e); %inicialización de la matriz de rigidez axil de cada elemento
    for ip = 1:IPKr
        chi = chi_pikr(ip); w = w_pikr(ip);
        Bef = (1/Je1)*[0 -1/2 0 1/2]; %B en flexión del elemento
        Bec = ((1/Je1)*[-1/2 0 1/2 0])-[0 (1-chi)/2 0 (1+chi)/2]; %B en cortadura del elemento
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        I_chi = NeL*I_nodalv(nodos_ev(g,:));
        A_chi = NeL*Av_c_nodal(nodos_ev(g,:));
        Kef = Kef + Bef'*E*I_chi*Bef*Je1*w;
        Kec = Kec + Bec'*G*A_chi*Bec*Je1*w;
        Kefc = Kef + Kec;
    end
    
    for ip = 1:IPK
        chi = chi_pik(ip); w = w_pik(ip);
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        Bea = (1/Je1)*[-1/2 1/2]; %elemento barra de dos nodos
        A_chi = NeL*A_nodalv(nodos_ev(g,:));
        Kea = Kea + Bea'*E*A_chi*Bea*Je1*w; 
    end
        Ke(gdl_flexion,gdl_flexion) = Kefc +  Ke(gdl_flexion,gdl_flexion); %matriz de rigidez completa del elemento con la contribución de flexión
        Ke(gdl_axil,gdl_axil) = Kea + Ke(gdl_axil,gdl_axil); %se añade la contribución axil
        Ke;
        K_quest(gdl_nodalv(g,:),gdl_nodalv(g,:))= Ke + K_quest(gdl_nodalv(g,:),gdl_nodalv(g,:));
end
K_quest;

%riostra
for g = 1:ne2
    Ke = zeros(nod_e*gdln); %inicialización de la matriz de rigidez en axil y flexión de cada elemento
    Kef = zeros(nod_e*2); %inicialización de la matriz de rigidez en flexión de cada elemento
    Kec = zeros(nod_e*2);
    Kefc = zeros(nod_e*2);
    Kea = zeros(nod_e); %inicialización de la matriz de rigidez axil de cada elemento
    for ip = 1:IPKr
        chi = chi_pikr(ip); w = w_pikr(ip);
        Bef = (1/Je2)*[0 -1/2 0 1/2]; %B en flexión del elemento
        Bec = ((1/Je2)*[-1/2 0 1/2 0])-[0 (1-chi)/2 0 (1+chi)/2]; %B en cortadura del elemento
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        I_chi = NeL*I_nodalr(nodos_er(g,:));
        A_chi = NeL*Ar_c_nodal(nodos_er(g,:));
        Kef = Kef + Bef'*E*I_chi*Bef*Je2*w;
        Kec = Kec + Bec'*G*A_chi*Bec*Je2*w;
        Kefc = Kef + Kec;
    end
    
    for ip = 1:IPK
        chi = chi_pik(ip); w = w_pik(ip);
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        Bea = (1/Je2)*[-1/2 1/2]; %elemento barra de dos nodos
        A_chi = NeL*A_nodalr(nodos_er(g,:));
        Kea = Kea + Bea'*E*A_chi*Bea*Je2*w; 
    end
        Ke(gdl_flexion,gdl_flexion) = Kefc +  Ke(gdl_flexion,gdl_flexion); %matriz de rigidez completa del elemento con la contribución de flexión
        Ke(gdl_axil,gdl_axil) = Kea + Ke(gdl_axil,gdl_axil); %se añade la contribución axil
        Ke;
        Ke_glob = Rr_'*Ke*Rr_;
        K_quest(gdl_nodalr(g,:),gdl_nodalr(g,:))= Ke_glob + K_quest(gdl_nodalr(g,:),gdl_nodalr(g,:));
end
K_quest

%% Matriz de coeficientes de restricción
p = 9; %número de restricciones totales del problema
r = zeros(p,1); %vector de restricciones del problema p = length(r)
C = zeros(p,size(K_quest,1)); %inicialización de la matriz de coeficientes de constraints

%SPCs (Single Point Constraint) - en los apoyos
C(1,1) = 1; C(2,2) = 1; C(3,17) = 1; C(4,19) = 1; C(5,20) = 1;

%SPCs (Single Point Constraint) - desplzamiento impuesto
C(6,14) = 1;

%MPCs (Multi Point Constraint) - unión no congruente
n = 2; %elemento en el que está la unión
x_chi = x_MPC-(n-1)*Le1; %x_chi es la coordenada x medida desde x1 del elemento en el que estamos
chi_ = 2*x_chi/Le1 - 1; %coordenada chi de la unión

n1 = (1-chi_)/2; n2 = (1+chi_)/2; %funciones lagrangianas para el desplazamiento axil en el elemento master

% h1 = -0.75*chi_+0.25*chi_^3+0.5; h2 = (-0.25*chi_-0.25*chi_^2+0.25*chi_^3+0.25)*Je1; %funciones hermíticas para el desplazamiento a flexión
% h3 = 0.75*chi_-0.25*chi_^3+0.5; h4 = (-0.25*chi_+0.25*chi_^2+0.25*chi_^3-0.25)*Je1; %elemento master
% h1_d = (-0.75+0.25*3*chi_^2)/Je1; h2_d = -0.25-0.25*2*chi_+0.25*3*chi_^2; %funciones hermíticas para el giro
% h3_d = (0.75-0.25*3*chi_^2)/Je1; h4_d = -0.25+0.25*2*chi_+0.25*3*chi_^2;

C(7,4) = n1; C(7,7) = n2; C(7,25) = -1;
C(8,5) = n1; C(8,8) = n2; C(8,26) = -1;
C(9,6) = n1; C(9,9) = n2; C(9,27) = -1;

r(6,1) = u5;

C_quest = C

%% Fuerza puntual equivalente al desplazamiento impuesto (N)
%MÉTODO DE PENALIZACIÓN
kpen = k*max(diag(K_quest)); %parámetro de penalización
k_ = kpen*eye(p);%matriz de penalización
u_pen = (K_quest+C_quest'*k_*C_quest)\(F+C_quest'*k_*r);
f_quest = K_quest*u_pen;
f_uload = f_quest(14) %coordenadas globales

% %MÉTODO DE LAGRANGE
% K_mL = [K_quest C_quest'; C_quest zeros(p)];%matriz de rigidez ampliada por los multiplicadores L
% F_mL = [F;r]; %vector de fuerzas ampliado
% U_mL = K_mL\F_mL; %vector de incógnitas: desplazamientos ampliados
% u_mL = U_mL(1:length(F)); %vector de desplzamientos nodales por m. de Lagrange
% lambda = U_mL(length(F)+1:end); %multiplicadores de Lagrange
% fR_mL = -lambda(1:5); %fuerzas de reacción por multiplicadores de Lagrange
% f_internas = -lambda(6:end); %fuerzas internas en la MPC por multiplicadores de Lagrange
% f_uload = f_internas(1);

%% Tensión de la fuerza axil sobre la riostra (MPa)
gdl_riostra = [nn1*gdln+1:1:nnt*gdln];
Rr_9 = [Rr zeros(3,6); zeros(3,3) Rr zeros(3,3); zeros(3,6) Rr]; 
f_r_glob = f_quest(gdl_riostra);
f_r_local = Rr_9*f_r_glob;
f_r = min(f_r_local);
sigma_axil_v2 = f_r/Ar*10^-6

%% Tensión cortante máxima sobre la viga principal, en valor absoluto (MPa)

z_max = a1/2;
for j = 1:ne1
    u_nodos_e = u_pen(gdl_nodalv(j, :));
    u_f_e = u_nodos_e(gdl_flexion); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en el elemento

    % deformaciones y tensiones en los puntos de integración
    ganma_ip_e = zeros(2,1); %deformación del elemento en los puntos de integración
    tau_ip_e = zeros(2,1);
    for ip = 1:IPKr
        chi = chi_pikr(ip);%datos de integración
        Bec = ((1/Je1)*[-1/2 0 1/2 0])-[0 (1-chi)/2 0 (1+chi)/2]; %B en cortadura del elemento
        ganma_ip_e(ip) = Bec*u_f_e;    
        tau_ip_e(ip) = ganma_ip_e(ip)*G;
        tau_ip(:,j) = tau_ip_e;
    end
    
end
tau_ip; %tensión cortante en los puntos de integración de cada elemento
tau_max_v1 = max(max(tau_ip))*10^-6 %tensión cortante máxima de la viga principal en megapascales

%% Tensión normal máxima de flexión sobre la viga principal, en valor absoluto (MPa)
z_max = a1/2;
for j = 1:ne1
    u_nodos_e = u_pen(gdl_nodalv(j, :));
    u_f_e = u_nodos_e(gdl_flexion); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en el elemento

    % deformaciones y tensiones en los puntos de integración
    epsilon_ip_e = zeros(2,1); %deformación del elemento en los puntos de integración
    sigma_ip_e = zeros(2,1);
    for ip = 1:IPKr
        chi = chi_pikr(ip);%datos de integración
        Bef = (1/Je1)*[0 -1/2 0 1/2]; %B en flexión del elemento
        epsilon_ip_e(ip) = z_max*Bef*u_f_e;    
        sigma_ip_e(ip) = epsilon_ip_e(ip)*E;
        sigma_ip(:,j) = sigma_ip_e;
    end
    
end
sigma_ip; %tensión en los puntos de integración de cada elemento

sigmaf_max_v1 = max(max(sigma_ip))*10^-6 %tensión máxima de la viga en megapascales


