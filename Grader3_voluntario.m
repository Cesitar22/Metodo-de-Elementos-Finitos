clear all;
close all;
clc;
%% Datos de entrada
E = 70e9; nu = 0.3; L = 1; a0 = L/25; b = 3*a0; t = a0/25; q0 = 2500; P = -1000;
ne = 10; % número de elementos viga hermíticos
nod_e = 2; % número de nodos de cada elemento (viga de dos nodos)
gdln = 3;%grados de libertad por nodo axil, cortante y flector
IPK = 3; %selección de la cantidad de puntos de integración por cada elemento, para la matriz de rigidez
IPF = 3; %selección de la cantidad de puntos de integración por cada elemento, para el vector de fuerzas nodales equivalentes
x_baricentro = 3/8*L; %baricentro de la carga distribuida semiparabólica
x_libre = L;
x_P = 5/10*L; %posición de la carga aplicada P
Le = L/ne; %longitud del elemento
n_P = round(x_P/Le); %nodo en el que está colocado la carga P

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

if IPF == 1;
    chi_pif = chi_pi1; w_pif = w_pi1;
elseif IPF == 2;
    chi_pif = chi_pi2; w_pif = w_pi2;
elseif IPF == 3;
    chi_pif = chi_pi3; w_pif = w_pi3;
else IPF == 4;
    chi_pif = chi_pi4; w_pif = w_pi4;
end

%% Definición de variables

Je = Le/2; iJe = 1/Je;%jacobiano e inverso del jacobiano
nn = ne +1; % número de nodos
gdlt = nn*gdln; %grados de libertad totales
gdl_T = [1:nn*gdln]; %gdl totales
gdl_R = [1, 2, 3]; %gdl restringidos
gdl_L = setdiff(gdl_T, gdl_R); %gdl libres

coord_n = zeros(nn, 3);
for i = 1:nn
    coord_n(i, 1) = (i-1)*Le;
end
coord_n; % matriz de coordenadas nodales: en filas nodos, en columna 1 coordenada x, en columna 2 coordenada y, en columna 3 coordenada z

fy_nodal = zeros(nn, 1);
for i = 1:nn
    x = coord_n(i,1);
    fy_nodal(i) = -q0/L^2*x^2 + q0;
end
fy_nodal; %valor de la fuerza cortante en los nodos

fx_nodal = zeros(nn, 1);
for i = 1:nn
    fx_nodal(i) = 0;
end
fx_nodal; %valor de la fuerza cortante en los nodos

nodos_e = zeros(ne, nod_e);
for i = 1:ne-1
    for j = 1:nod_e
        nodos_e(1,j) = j;
        nodos_e(i+1,j) = nodos_e(1,j)+i;
    end
end
nodos_e; %nodos de cada elemento (una fila es un elemento)

z_x = zeros(nn,1);
x_n = [0:Le:L]; %coordenadas de los nodos
for i = 1:nn
    x_i = x_n(i);
    a_x = (-2/3*a0/L)*x_i + a0;
    z_x(i) = a_x/2;
end
z_x; %vector de 'z' en los nodos

Z_x = zeros(ne,nod_e);
for i = 1:ne
    Z_x (i,[1:nod_e]) = z_x(nodos_e(i,[1:nod_e]));
end
Z_x; %matriz de valores 'z' en los nodos, un elemento por fila

x_medio = [Le/2:Le:L-Le/2]'; %coordenadas de los puntos medios de los elementos

x__1 = [0:Le:L-Le]; %posición del primer nodos de cada elemento
x_p1 = zeros(ne,1); %posición del primer punto de integración de cada elemento
x_p2 = zeros(ne,1); %posición del segundo punto de integración de cada elemento
for i = 1:ne
    y_p1 = (-1/sqrt(3) + 1)*Le/2;
    y_p2 = (1/sqrt(3) + 1)*Le/2;
    x_p1(i) = x__1(i) + y_p1;
    x_p2(i) = x__1(i) + y_p2;
end
x_int = [x_p1, x_p2]; %coordenadas de los puntos de integración de todos los elementos

z_x_int = zeros(ne,2);
for i = 1:ne
    for j =1:2
    x_i = x_int(i,j);
    a_x = (-2/3*a0/L)*x_i + a0;
    z_x_int(i,j) = a_x/2;
    end
end
z_x_int; %vector de 'z' en los puntos de integración, una fila por cada elemento, siendo las columnas el punto de integración 1 y el 2

gdl_nodal = zeros(ne,gdln*2); 
for i = 1:ne-1
   for j = 1:gdln*2
           gdl_nodal(1,j) = j; 
           gdl_nodal(i+1,j) = gdl_nodal(1,j)+3*i;
    end
end
gdl_nodal; %grados de libertad que le corresponden a cada a ELEMENTO

gdl_nod = zeros(nn,gdln);
for i = 1:nn-1
   for j = 1:gdln
           gdl_nod(1,j) = j; 
           gdl_nod(i+1,j) = gdl_nod(1,j)+3*i;
    end
end
gdl_nod; %grados de libertad que le corresponden a cada a NODO


I_nodal = zeros(nn, 1);
for i = 1:nn
    x = coord_n(i,1);
    a_x = (-2/3*a0/L)*x + a0;
    b_x = b;
    t_x = t;
    I_nodal(i) = 1/12*a_x^3*b_x-1/12*(a_x-2*t_x)^3*(b_x-2*t_x); %calculada como el momento de inercia del rectángulo exterior menos el rectángulo interior
end
I_nodal; %valor de la fuerza cortante en los nodos

A_nodal = zeros(nn, 1);
for i = 1:nn
    x = coord_n(i,1);
    a_x = (-2/3*a0/L)*x + a0;
    b_x = b;
    t_x = t;
    A_nodal(i,1) = a_x*b_x-(a_x-2*t_x)*(b_x-2*t_x);
end
A_nodal; %valor del área en los nodos


%% Cálculo de K_global
K_quest = zeros(gdln*nn);
gdl_flexion = [2, 3, 5, 6];
gdl_axil = [1 4];
for g = 1:ne
    Ke = zeros(nod_e*gdln); %inicialización de la matriz de rigidez en axil y flexión de cada elemento
    Kef = zeros(nod_e*2); %inicialización de la matriz de rigidez en flexión de cada elemento
    Kea = zeros(nod_e); %inicialización de la matriz de rigidez axil de cada elemento
    for ip=1:IPK
        chi = chi_pik(ip); w = w_pik(ip);
        k1 = (chi*3/2); k2 = (chi*3-1)/2; k3 = -(chi*3/2); k4 = (chi*3+1)/2; %curvaturas de la viga
        Bef = iJe^2*[k1 Je*k2 k3 Je*k4]; %elemento viga de 2 nodos
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        I_chi = NeL*I_nodal(nodos_e(g,:)); 
        Kef = Kef + Bef'*E*I_chi*Bef*Je*w;
        
        Bea = (1/Je)*[-1/2 1/2]; %elemento barra de dos nodos
        A_chi = NeL*A_nodal(nodos_e(g,:));
        Kea = Kea + Bea'*E*A_chi*Bea*Je*w;
    end
        Ke(gdl_flexion,gdl_flexion) = Kef +  Ke(gdl_flexion,gdl_flexion); %matriz de rigidez completa del elemento con la contribución de flexión
        Ke(gdl_axil,gdl_axil) = Kea + Ke(gdl_axil,gdl_axil); %se añade la contribución axil
        Ke;
        K_quest(gdl_nodal(g,:),gdl_nodal(g,:))= Ke + K_quest(gdl_nodal(g,:),gdl_nodal(g,:));
end
K_quest



%% Cálculo de vector de fuerzas de flexión nodales equivalentes
F_quest = zeros(gdln*nn,1);
F_P = zeros(gdln*nn,1);
F_P(gdln*(n_P-1)+2) = P;
for g = 1:ne
    Fe = zeros(nod_e*gdln,1); %inicialización del vector de fuerzas nodales total de cada elemento
    Fef = zeros(nod_e*2,1); %inicialización del vector de fuerzas nodales en flexión de cada elemento
    Fea = zeros(nod_e,1); %inicialización del vector de fuerzas nodales axil de cada elemento
    for ip=1:IPF
        chi = chi_pif(ip); w = w_pif(ip);
        nh1 = -0.75*chi+0.25*chi^3+0.5; nh2 = -0.25*chi-0.25*chi^2+0.25*chi^3+0.25;
        nh3 = 0.75*chi-0.25*chi^3+0.5; nh4 = -0.25*chi+0.25*chi^2+0.25*chi^3-0.25;
        Neh = [nh1 Je*nh2 nh3 Je*nh4];
        fynod = fy_nodal(nodos_e(g,:)); %valor en los nodos de la fuerza cortante externa distribuida
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        fy_chi = NeL*fynod; %fuerza externa virgulilla (en formulación isoparaméttrica)
        Fef = Fef + Neh'*fy_chi*Je*w;
        
        f_chi = NeL*fx_nodal(nodos_e(g,:)); %fuerza externa virgulilla (en formulación isoparaméttrica)
        Fea = Fea + NeL'*f_chi*Je*w;
    end
        Fe(gdl_flexion,1) = Fef +  Fe(gdl_flexion,1); %matriz de rigidez completa del elemento con la contribución de flexión
        Fe(gdl_axil,1) = Fea + Fe(gdl_axil,1); %se añade la contribución axil
        Fe;
        F_quest(gdl_nodal(g,:),1)= Fe + F_quest(gdl_nodal(g,:),1);
end
F_quest = F_quest + F_P;

%% Cálculo de desplazamientos en el extremo libre y en el baricentro

KLL = K_quest(gdl_L, gdl_L); %matriz de rigidez reducida
F_L = F_quest(gdl_L);
uL_quest = KLL\F_L; %vector de desplazamientos en nodos libres
U = zeros(length(gdl_T),1);
U(gdl_L) = uL_quest; %vector de desplazamientos completo
U;

u1_quest = U(gdl_nod(end,:),1); %desplazamientos obtenidos en el extremo libre
u1_quest = round(u1_quest, 4)

%% Cálculo de tensiones; solución fuerte
x_n1 = [0:Le:L-Le];
x_n2 = [Le:Le:L];
sigma_x = zeros(length(x_n1), nod_e);
for n = 1:length(x_n1)
    x_i = x_n1(n);
    a_x = (-2/3*a0/L)*x_i + a0;
    b_x = b;
    t_x = t;
    
    z_max = a_x/2;
    I_y = 1/12*a_x^3*b_x-1/12*(a_x-2*t_x)^3*(b_x-2*t_x);
    My_x = q0*L^2*(-0.25 + (2/3)*(x_i/L) - 0.5*(x_i/L)^2 + (1/12)*(x_i/L)^4);
    sigma_x(n,1) = abs(-My_x/I_y*z_max);
end
for n = 1:length(x_n2)
    x_i = x_n2(n);
    a_x = (-2/3*a0/L)*x_i + a0;
    b_x = b;
    t_x = t;
    
    z_max = a_x/2;
    I_y = 1/12*a_x^3*b_x-1/12*(a_x-2*t_x)^3*(b_x-2*t_x);
    My_x = q0*L^2*(-0.25 + (2/3)*(x_i/L) - 0.5*(x_i/L)^2 + (1/12)*(x_i/L)^4);
    sigma_x(n,2) = abs(-My_x/I_y*z_max);
end
sigma_x; %tensión en los nodos de cada elemento
sigma1_fuerte = max(max(sigma_x)); %tensión máxima de la viga, y su posición

%% Cálculo de tensiones en los puntos de integración y los nodos; solución débil
for j = 1:ne
    u_nodos_e = U(gdl_nodal(j, :));
    u_f_e = u_nodos_e(gdl_flexion); % define las componentes o grados de libertad que forman parte del vector de desplazamientos en el elemento

    % deformaciones y tensiones en los puntos de integración
    epsilon_ip_e = zeros(2,1); %deformación del elemento en los puntos de integración
    sigma_ip_e = zeros(2,1);
    for ip = 1:2
        chi = chi_pi2(ip);%datos de integración
        nh1_2d = 3/2*chi; nh2_2d = -1/2+3/2*chi;
        nh3_2d = -3/2*chi; nh4_2d = 1/2+3/2*chi;
        Neh_2d = iJe^2*[nh1_2d Je*nh2_2d nh3_2d Je*nh4_2d]; %Bf
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        z = NeL*(Z_x(j,:))';
        epsilon_ip_e(ip) = z*Neh_2d*u_f_e;    
        sigma_ip_e(ip) = epsilon_ip_e(ip)*E;
        sigma_ip(:,j) = sigma_ip_e;
    end
    
end
sigma_ip; %tensión en los puntos de integración de cada elemento
for j = 1:ne
    for i = 1:2
        chi_p = [-sqrt(3) sqrt(3)];
        chi = chi_p(i);
        n1L = (1-chi)/2; n2L = (1+chi)/2;
        NeL = [n1L n2L]; %funciones de interpolación isoparamétrica para dos nodos: lagrangianas
        sigma_nod(i,j) = NeL*sigma_ip(:,j);
    end
end
sigma_nod; %tensión en los nodos, calculada con la interpolación de Lagrange
sigma1_quest = max(max(sigma_nod))/(10^6) %tensión en megapascales

% %% Representación gráfica de tensiones normales y cortantes de flexión
% % Esta representación se puede hacer en bucle igual que antes la integración numérica y el ensamblaje
% % Calcula las tensiones en cada elemento
% Mze_i = zeros(2,ne) %matriz de almacen de momentos en los nodos, en columnas por elementos y en filas por puntos de integración
% sigmaf_i = zeros(2,ne) %matriz de almacen de sigma en los nodos, en columnas por elementos y en filas por puntos de integración
% for e = 1:ne 
%     index = conectividad_e(e,:); %indice de posición de los nodos en global
%     x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
%     %...
%        
%     %calcula la tension en los nodos de cada elemento y elemento a elemento
%     Mze_i(:,e)
%     sigmaf_i(:,e) =
%     
%     xi_ip = % coordendas de los puntos de integración en coordenadas globales: usar la transformación isoparamétrica en geometría x~ = N*x^
%     figure(1)   
%     graf1 = plot(xi_ip,Mze_i,'bs-','LineWidth',1); % xi_ip son las coordenadas de los puntos de integración en coordenadas globales,
%     axis([0 1 0 max(Mze_i)*1.02]); % ajuste de la escala de ejes de representación 
%     title("Momento Flecor Nm",'Fontsize',12)
%     xlabel("Nm"); ylabel("m");
%     hold on
%     
%     figure(2)
%     graf2 = plot(xi_ip,sigmaf_i,'bs-','LineWidth',1); % xi_ip son las coordenadas de los puntos de integración en coordenadas globales, y sigmafi es la tensión de flexión en MPa
%     axis([0 1 0 max(sigmaf_i)*1.02]); % ajuste de la escala de ejes de representación con sigmafi en MPa
%     title("Tensión normal MPa",'Fontsize',12)
%     xlabel("MPa"); ylabel("m");
%     hold on
%     
% end

%% Representa desplazamientos
GDL = gdln*nn;
Ux = U(1:3:(GDL));
Uy = U(2:3:(GDL)); Uymax = max(abs(Uy)); %flecha maxima
Uxy = U(3:3:(GDL));
escala = 0.15*L/Uymax;
coordx_n = coord_n(:,1);
coordy_n = coord_n(:,2);
graf3 = plot(coordx_n,coordy_n,'bo-',coordx_n+Ux,coordy_n+Uy,'rs--','LineWidth',2)
axis([0 (L+max(Ux)*10) min(Uy)*2 max(Uy)*2])
title("Desplazamientos nodales",'Fontsize',12)
xlabel("coordenada nodal X"); ylabel("coordenada nodal Y");

% %% Opcional: si quieres representar las tensiones nodales (en otra figura)
% for e = 1:ne
%     
% sigmaf_n = zeros(e;2) %matriz de almacen de sigma en los nodos, en filas por elementos y en columnas por nodos
% Mze_n = zeros(e;2) %matriz de almacen de sigma en los nodos, en filas por elementos y en columnas por nodos
% 
% %calcula la tension en los nodos de cada elemento y elemento a elemento
% sigmaf_n(e,:) =
% Mze_n(e,:) =
% 
% figure(3)
% plot(x_n,Mze_n(e,:),'ro--','LineWidth',1); % x_n son las coordenadas de los nodos, y tau_n es la tensión en MPa
% hold on
% 
% figure(4)
% plot(x_n,sigmaf_n(e,:),'ro--','LineWidth',1);   % x_n son las coordenadas de los nodos, y sigmaf_n es la tensión en MPa
% hold on
% 
% end