close all
clear all
clc

% Inicio del código
L = 1; q0 = 15000; E = 1000; I = 10; A = 100;
%longitud de los elementos
Le1 = L; Le2 = L; Le3 = sqrt(L^2 + (L/2)^2);
Le4 = sqrt(L^2 + (L/2)^2); Le5 = L/2; Le6 = L/2; 
Le7 = sqrt(L^2 + (L/2)^2);Le8 = sqrt(L^2 + (L/2)^2);
%orientación de los elementos
alfa1 = pi/2; alfa2 = pi/2; alfa3 = -atan((L/2)/L);
alfa4 = -atan((L/2)/L); alfa5 = -pi/2; alfa6 = -pi/2;
alfa7 = atan((L/2)/L); alfa8 = pi + atan((L/2)/L);

zerosR = zeros(3);

%elemento 1
Ke1 = E*I*[A*E/(Le1*E*I) 0 0 -A*E/(Le1*E*I) 0 0; 0 12/(Le1^3) 6/(Le1^2) 0 -12/(Le1^3) 6/(Le1^2); 0 6/(Le1^2) 4/Le1 0 -6/(Le1^2) 2/Le1; -A*E/(Le1*E*I) 0 0 A*E/(Le1*E*I) 0 0; 0 -12/(Le1^3) -6/(Le1^2) 0 12/(Le1^3) -6/(Le1^2); 0 6/(Le1^2) 2/Le1 0 -6/(Le1^2) 4/Le1]; %matriz de rigidez local del elemento
Re1 = [cos(alfa1) sin(alfa1) 0; -sin(alfa1) cos(alfa1) 0; 0 0 1];
Re1_ = [Re1 zerosR; zerosR Re1];
Ke1_ = Re1_'*Ke1*Re1_; %matriz de rigidez global del elemento
%elemento 2
Ke2 = E*I*[A*E/(Le2*E*I) 0 0 -A*E/(Le2*E*I) 0 0; 0 12/(Le2^3) 6/(Le2^2) 0 -12/(Le2^3) 6/(Le2^2); 0 6/(Le2^2) 4/Le2 0 -6/(Le2^2) 2/Le2; -A*E/(Le2*E*I) 0 0 A*E/(Le2*E*I) 0 0; 0 -12/(Le2^3) -6/(Le2^2) 0 12/(Le2^3) -6/(Le2^2); 0 6/(Le2^2) 2/Le2 0 -6/(Le2^2) 4/Le2]; %matriz de rigidez local del elemento
Re2 = [cos(alfa2) sin(alfa2) 0; -sin(alfa2) cos(alfa2) 0; 0 0 1];
Re2_ = [Re2 zerosR; zerosR Re2];
Ke2_ = Re2_'*Ke2*Re2_; %matriz de rigidez global del elemento
%elemento 3
Ke3 = E*I*[A*E/(Le3*E*I) 0 0 -A*E/(Le3*E*I) 0 0; 0 12/(Le3^3) 6/(Le3^2) 0 -12/(Le3^3) 6/(Le3^2); 0 6/(Le3^2) 4/Le3 0 -6/(Le3^2) 2/Le3; -A*E/(Le3*E*I) 0 0 A*E/(Le3*E*I) 0 0; 0 -12/(Le3^3) -6/(Le3^2) 0 12/(Le3^3) -6/(Le3^2); 0 6/(Le3^2) 2/Le3 0 -6/(Le3^2) 4/Le3]; %matriz de rigidez local del elemento
Re3 = [cos(alfa3) sin(alfa3) 0; -sin(alfa3) cos(alfa3) 0; 0 0 1];
Re3_ = [Re3 zerosR; zerosR Re3];
Ke3_ = Re3_'*Ke3*Re3_; %matriz de rigidez global del elemento
%elemento 4
Ke4 = E*I*[A*E/(Le4*E*I) 0 0 -A*E/(Le4*E*I) 0 0; 0 12/(Le4^3) 6/(Le4^2) 0 -12/(Le4^3) 6/(Le4^2); 0 6/(Le4^2) 4/Le4 0 -6/(Le4^2) 2/Le4; -A*E/(Le4*E*I) 0 0 A*E/(Le4*E*I) 0 0; 0 -12/(Le4^3) -6/(Le4^2) 0 12/(Le4^3) -6/(Le4^2); 0 6/(Le4^2) 2/Le4 0 -6/(Le4^2) 4/Le4]; %matriz de rigidez local del elemento
Re4 = [cos(alfa4) sin(alfa4) 0; -sin(alfa4) cos(alfa4) 0; 0 0 1];
Re4_ = [Re4 zerosR; zerosR Re4];
Ke4_ = Re4_'*Ke4*Re4_; %matriz de rigidez global del elemento
%elemento 5
Ke5 = E*I*[A*E/(Le5*E*I) 0 0 -A*E/(Le5*E*I) 0 0; 0 12/(Le5^3) 6/(Le5^2) 0 -12/(Le5^3) 6/(Le5^2); 0 6/(Le5^2) 4/Le5 0 -6/(Le5^2) 2/Le5; -A*E/(Le5*E*I) 0 0 A*E/(Le5*E*I) 0 0; 0 -12/(Le5^3) -6/(Le5^2) 0 12/(Le5^3) -6/(Le5^2); 0 6/(Le5^2) 2/Le5 0 -6/(Le5^2) 4/Le5]; %matriz de rigidez local del elemento
Re5 = [cos(alfa5) sin(alfa5) 0; -sin(alfa5) cos(alfa5) 0; 0 0 1];
Re5_ = [Re5 zerosR; zerosR Re5];
Ke5_ = Re5_'*Ke5*Re5_; %matriz de rigidez global del elemento
%elemento 6
Ke6 = E*I*[A*E/(Le6*E*I) 0 0 -A*E/(Le6*E*I) 0 0; 0 12/(Le6^3) 6/(Le6^2) 0 -12/(Le6^3) 6/(Le6^2); 0 6/(Le6^2) 4/Le6 0 -6/(Le6^2) 2/Le6; -A*E/(Le6*E*I) 0 0 A*E/(Le6*E*I) 0 0; 0 -12/(Le6^3) -6/(Le6^2) 0 12/(Le6^3) -6/(Le6^2); 0 6/(Le6^2) 2/Le6 0 -6/(Le6^2) 4/Le6]; %matriz de rigidez local del elemento
Re6 = [cos(alfa6) sin(alfa6) 0; -sin(alfa6) cos(alfa6) 0; 0 0 1];
Re6_ = [Re6 zerosR; zerosR Re6];
Ke6_ = Re6_'*Ke6*Re6_; %matriz de rigidez global del elemento
%elemento 7
Ke7 = E*I*[A*E/(Le7*E*I) 0 0 -A*E/(Le7*E*I) 0 0; 0 12/(Le7^3) 6/(Le7^2) 0 -12/(Le7^3) 6/(Le7^2); 0 6/(Le7^2) 4/Le7 0 -6/(Le7^2) 2/Le7; -A*E/(Le7*E*I) 0 0 A*E/(Le7*E*I) 0 0; 0 -12/(Le7^3) -6/(Le7^2) 0 12/(Le7^3) -6/(Le7^2); 0 6/(Le7^2) 2/Le7 0 -6/(Le7^2) 4/Le7]; %matriz de rigidez local del elemento
Re7 = [cos(alfa7) sin(alfa7) 0; -sin(alfa7) cos(alfa7) 0; 0 0 1];
Re7_ = [Re7 zerosR; zerosR Re7];
Ke7_ = Re7_'*Ke7*Re7_; %matriz de rigidez global del elemento
%elemento 8
Ke8 = E*I*[A*E/(Le8*E*I) 0 0 -A*E/(Le8*E*I) 0 0; 0 12/(Le8^3) 6/(Le8^2) 0 -12/(Le8^3) 6/(Le8^2); 0 6/(Le8^2) 4/Le8 0 -6/(Le8^2) 2/Le8; -A*E/(Le8*E*I) 0 0 A*E/(Le8*E*I) 0 0; 0 -12/(Le8^3) -6/(Le8^2) 0 12/(Le8^3) -6/(Le8^2); 0 6/(Le8^2) 2/Le8 0 -6/(Le8^2) 4/Le8]; %matriz de rigidez local del elemento
Re8 = [cos(alfa8) sin(alfa8) 0; -sin(alfa8) cos(alfa8) 0; 0 0 1];
Re8_ = [Re8 zerosR; zerosR Re8];
Ke8_ = Re8_'*Ke8*Re8_; %matriz de rigidez global del elemento

%ensamblaje
K_quest = zeros(24);

gdl_e1 = [1 2 3 4 5 6];
K_quest(gdl_e1, gdl_e1) = Ke1_ + K_quest(gdl_e1, gdl_e1);
gdl_e2 = [4 5 6 7 8 9];
K_quest(gdl_e2, gdl_e2) = Ke2_ + K_quest(gdl_e2, gdl_e2);
gdl_e3 = [7 8 9 10 11 12];
K_quest(gdl_e3, gdl_e3) = Ke3_ + K_quest(gdl_e3, gdl_e3);
gdl_e4 = [10 11 12 13 14 15];
K_quest(gdl_e4, gdl_e4) = Ke4_ + K_quest(gdl_e4, gdl_e4);
gdl_e5 = [13 14 15 16 17 18];
K_quest(gdl_e5, gdl_e5) = Ke5_ + K_quest(gdl_e5, gdl_e5);
gdl_e6 = [16 17 18 19 20 21];
K_quest(gdl_e6, gdl_e6) = Ke6_ + K_quest(gdl_e6, gdl_e6);
gdl_e7 = [1 2 3 22 23 24];
K_quest(gdl_e7, gdl_e7) = Ke7_ + K_quest(gdl_e7, gdl_e7);
gdl_e8 = [13 14 15 22 23 24];
K_quest(gdl_e8, gdl_e8) = Ke8_ + K_quest(gdl_e8, gdl_e8);

K_quest = round(K_quest) %matriz de rigidez global

%condiciones de contorno
gdl_L = [[4:1:18] [22 23 24]];%grados de libertad libres
gdl_R = [1 2 3 19 20 21]; %grados de libertad restringidos
KLL = K_quest(gdl_L, gdl_L); %submatriz de rigidez de grados libres
KRR = K_quest(gdl_R, gdl_R); %submatriz de rigidez de grados restringidos
KLR = K_quest(gdl_L, gdl_R); %términos cruzados
KRL = KLR'; %términos cruzados

% vector de fuerzas nodales equivalentes en coord. locales
f1 = [0 0 0]';
f2 = [0 0 0]';
f3 = [0 -3/20*q0/2*Le3 -q0/2*(Le3^2)/30]';%en coordenadas locales del elemento 3-4
f4 = [0 -7/20*q0/2*Le3-13/40*q0*Le4 q0/2*(Le3^2)/20-7/120*q0*Le4^2]';%en cordenadas locales del elemento 3-4
f5 = [0 -17/40*Le4*q0 q0*(Le4^2)/15]'; %en coordenadas locales del elemento 3-4
f6 = [0 0 0]';
f7 = [0 0 0]';
f8 = [0 0 0]';

% vector de fuerzas nodales equivalentes en coord. globales
f1_glob = Re1'*f1;
f2_glob = Re2'*f2;
f3_glob = Re3'*f3;
f4_glob = Re4'*f4;
f5_glob = Re4'*f5;
f6_glob = Re6'*f6;
f7_glob = Re7'*f7;
f8_glob = Re8'*f8;

feq_quest = [f1_glob; f2_glob; f3_glob; f4_glob; f5_glob; f6_glob; f7_glob; f8_glob] %vector de fuerzas nodales equivalentes de la estructura en ejes globales


% vector de desplazamientos nodales en los nodos centrales
fL = feq_quest(gdl_L);
uL = KLL\fL;
gdl_centrales = [1 2 3 7 8 9 13 14 15 16 17 18];
u_quest = uL(gdl_centrales)

% fuerza de reacción en apoyo izquierdo
fR = KRL*uL;
fr_quest = fR([1 2 3])

%vector de fuerzas internas locales del elemento 7
u_18 = [[0; 0; 0]; uL([16 17 18])]; %desplazamiento de los nodos 1 y 8 en coordenadas globales
u_18_ = Re7_*u_18; %desplazamiento de los nodos 1 y 8 en coordenadas locales
fele_quest = Ke7*u_18_

% Fin del código