%% Datos

v_infinito = 1
alpha = 10.0*pi/180.0 % radianes
rho = 1.29  %densidad en kg/m^3 del aire

%% Geometria

syms x
%Dijo en clase que lo mejor seria si hubiese mas nodos por los extremos, y
%dividir el perfil en angulos iguales, si nos resulta muy jodido siempre lo
%podemos hacer equispaciado

%Escribimos la ecuacion del perfil 
%Media aritmetica de nuestos DNI: 66325243. Por tanto nuestro perfil es
%NACA 5243

%Codio para dividir el perfil en angulos iguales
n_alpha = 10 % N + 1 Nodos
alpha1 = linspace(0, 2*pi, n_alpha); %tiene que repetirse un nodo en x=1
x_nodos = zeros(1,n_alpha); %Componente x de los nodos
for i=1:n_alpha
    x_nodos(i) = 0.5*(1 + cos(alpha1(i))); %Sacamos la componente x de los nodos
end

%Ahora escribimos la ecuacion del perfil
%Media aritmetica de nuestos DNI: 66325243. Por tanto nuestro perfil es
%NACA 5243

fmax = (1 + 2)*0.01; %Lo multiplico por 0.01 para estaba en porcentaje
x_fmax = (30 + 10*0)*0.01; %Lo multiplico por 0.1 porq
tmax = (7 + 4)*0.01; %Lo multiplico por 0.01 para estaba en porcentaje

z_curv1 = (2*x_fmax*x - x^2)*fmax/x_fmax^2; % (0<x<x_fmax) 
z_curv2 = ((1 - 2*x_fmax) + 2*x_fmax*x - x^2)*fmax/(1 - x_fmax)^2; % (x_fmax<x<1)

z_esp = 5*tmax*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4); %Hay que sumarla o restarla si es el extrados o el intrados respectivamente

z_extra1 = z_curv1 + z_esp ;% (0<x<x_fmax) 
z_extra2 = z_curv2 + z_esp; % (x_fmax<x<1)
z_intra1 = z_curv1 - z_esp; % (0<x<x_fmax) 
z_intra2 = z_curv2 - z_esp ;% (x_fmax<x<1)

%Sacamos la componente z de los nodos
%Pero primero tenemos que sacar el nº de nodos que estan en el extrados que
%es el mismo nº de los que estan en el intrados, dentro del contador_extra
%se cuenta el primer nodo cuando z=0 y x=1

if mod(n_alpha,2) == 0 % numero de nodos par, numero de paneles impar
    contador_extra = n_alpha/2; %numero de puntos en el extrados
else
    contador_extra = (n_alpha - 1)/2;  % En el caso de que el numero de nodos sea impar, y nº de paneles par
end 
z_nodos = zeros(1,n_alpha) ;
z_nodos(1) = 0;
z_nodos(n_alpha) = 0;
for i=2:n_alpha - 1
    if x_nodos(i)==0 %Solo en el caso de que n_alpha sea impar
        z_nodos(i) = 0;
    else
        if x_nodos(i)< x_fmax %1
            if i <= contador_extra %extrados
                z_nodos(i) = subs(z_extra1, x, x_nodos(i));
            else %intrados
                z_nodos(i) = subs(z_intra1, x, x_nodos(i));
            end
        else %2 
            if i <= contador_extra %extrados
                z_nodos(i) = subs(z_extra2, x, x_nodos(i));
            else %intrados
                z_nodos(i) = subs(z_intra2, x, x_nodos(i));
            end
        end
    end
end

%% Código

N = n_alpha - 1; %Numero de paneles
theta_1 = zeros(1, N); % vector fila
long_panel = zeros(1,N);
for i=1: N  % Hay N paneles
      theta_1(i) = atan2((z_nodos(i + 1) - z_nodos(i)),(x_nodos(i + 1) - x_nodos(i))) ;  % (0<theta<2*pi) comprobarlo y sino usar atan2 y sumarle 2*pi cuando sea negativo para que salga dentro de ese intervalo
      long_panel(i) = sqrt((x_nodos(i + 1) - x_nodos(i))^2 + (z_nodos(i + 1) - z_nodos(i))^2);
end 
theta = zeros(1, N);
for i=1:N
    if theta_1(i)>0
        theta(i) = theta_1(i);
    else
        theta(i) = 2*pi + theta_1(i);
    end
end
teta1 = theta_1(1)
teta4= theta_1(4)

long_panel1 = long_panel
theta1 = theta*180.0/pi
p = pi


x_c_panel = zeros(1,N); %coordenada x absoluta del centro del panel 
z_c_panel = zeros(1,N); %coordenada z absoluta del centro del panel
for i=1:N
    x_c_panel(i) = 0.5*(x_nodos(i+1) + x_nodos(i));
    z_c_panel(i) = 0.5*(z_nodos(i+1) + z_nodos(i));
end
x_c_panel1 = x_c_panel
z_c_panel1 = z_c_panel

x_c_panel_prima = zeros(N); %coordenada x relativa del centro del panel 
z_c_panel_prima = zeros(N); %coordenada z relativa del centro del panel 
beta = zeros(N);
for i=1:N %recorre el centro de todos los paneles, N
    for j=1:N  %Recorre todos los nodos pero no hace falta que llegue hasta el ultimo
        x_c_panel_prima(i,j) = (x_c_panel(i) - x_nodos(j))*cos(theta(j)) + (z_c_panel(i) - z_nodos(j))*sin(theta(j));
        z_c_panel_prima(i,j) = (z_c_panel(i) - z_nodos(j))*cos(theta(j)) - (x_c_panel(i) - x_nodos(j))*sin(theta(j));

        beta(i,j) = atan2(z_c_panel_prima(i,j),(x_c_panel_prima(i,j) - long_panel(j))) - atan2(z_c_panel_prima(i,j),x_c_panel_prima(i,j));
    end
end
Beta1 = beta %esta bien y por tanto tambien xc panel prima
r = zeros(N, N + 1);  
for i=1:N %recorre el centro de todos los paneles, N
    for j=1:N + 1 %Recorre todos los nodos N + 1
        r(i,j) = sqrt((x_c_panel(i)-x_nodos(j))^2 + (z_c_panel(i)-z_nodos(j))^2);
    end
end 
a5 = zeros(N);
a3 = zeros(N);
a6 = zeros(N);
a4 = zeros(N);
a7 = zeros(N);
a8 = zeros(N);
a9 = zeros(N);
a10 = zeros(N);
for i=1:N %Recorre todos los paneles, N
    for j=1:N
        a5(i,j) = (beta(i,j) - ((log(r(i,j + 1)/r(i,j))*z_c_panel_prima(i,j) + beta(i,j)*x_c_panel_prima(i,j))/long_panel(j)))/(2*pi);
        a3(i,j) = (log(r(i,j + 1)/r(i,j))*z_c_panel_prima(i,j) + beta(i,j)*x_c_panel_prima(i,j))/(long_panel(j)*2*pi);
        a6(i,j) = (log(r(i,j + 1)/r(i,j)) - ((log(r(i,j + 1)/r(i,j))*x_c_panel_prima(i,j) + long_panel(j) - beta(i,j)*z_c_panel_prima(i,j))/long_panel(j)))/(2*pi);
        a4(i,j) = (log(r(i,j + 1)/r(i,j))*x_c_panel_prima(i,j) + long_panel(j) - beta(i,j)*z_c_panel_prima(i,j))/(long_panel(j)*2*pi);
        a7(i,j) = a5(i,j)*cos(theta(j)) - a6(i,j)*sin(theta(j));
        a8(i,j) = a3(i,j)*cos(theta(j)) - a4(i,j)*sin(theta(j));
        a9(i,j) = a5(i,j)*sin(theta(j)) + a6(i,j)*cos(theta(j));
        a10(i,j) = a3(i,j)*sin(theta(j)) + a4(i,j)*cos(theta(j));
    end
end

B = zeros(N, N + 1); % N filas recorrida por i, N + 1 columnas recorrida por j
C = zeros(N, N + 1); % N filas recorrida por i, N + 1 columnas recorrida por j
A = zeros(N + 1);
b = zeros(1,N + 1);

for i=1:N 
    B(i,1) = a7(i,1); %matriz(fila, columna)
    B(i,N + 1) = a8(i,N);  %esta bien
    C(i,1) = a9(i,1); 
    C(i,N + 1) = a10(i,N); %esta bien
    b(i) = -v_infinito*sin(alpha -  theta(i));
    for j=2: N 
        B(i,j) = a8(i,j - 1) + a7(i,j);
        C(i,j) = a10(i,j - 1) + a9(i,j);
    end
end

for i=1: N 
    for j=1: N + 1
         A(i,j) = -B(i,j)*sin(theta(i)) + C(i,j)*cos(theta(i));
    end
end

b(N + 1) = 0;
A(N + 1,1) = 1;
A(N + 1,N + 1) = 1;
for j=2:N 
    A(N + 1,j) = 0;
end
B1 = B 
C1 = C 
A1 = A

gamma = A\b' % Resolvemos el sistema de ecuaciones lineal 

v_tan_panel = zeros(1,N); %velocidad tangencial de cada panel
for i=1: N
    v_tan_inducida = 0;
    for j =1: N + 1
        v_tan_inducida = v_tan_inducida +gamma(j)*(B(i,j)*cos(theta(i)) + C(i,j)*sin(theta(i))); %velocidad tangencial de cada panel inducida por los torbellinos
    end
    v_tan_panel(i) = v_infinito*cos(theta(i) - alpha) + v_tan_inducida;
end

Cp = zeros(1,N); %Coeficiente de presiones de cada panel
for i=1: N 
    Cp(i) = 1 - (v_tan_panel(i)/v_infinito)^2;
end

CP1 = Cp

Cl = zeros(1, N); %Coeficiente de sustentacion en cada panel
for i=1: N
    Cl(i) = (gamma(i) + gamma(i + 1))/v_infinito;
end
Cl_paneles = Cl 
CL_total = 0; %Coeficiente de sustentacion total 
for i=1:N 
    CL_total = CL_total +  Cl(i)*long_panel(i); 
end
CL = CL_total

l = zeros(1, N); %sustentacion de cada panel
for i=1: N
    l(i) = rho*v_infinito*long_panel(i)*(gamma(i) + gamma(i + 1))/2;
end

m_ab = zeros(1, N); % momento de la sustentacion de cada panel desde el borde de ataque
d = zeros(1, N); %Brazo del momento de la sustentacion
delta = zeros(1, N); %angulo que forma el brazo con el eje x
for i=1: N 
    d(i) = sqrt(x_c_panel(i)^2.0 + z_c_panel(i)^2.0);
    delta(i) = atan(z_c_panel(i)/x_c_panel(i));
    m_ab(i) = -l(i)*d(i)*cos(delta(i) - alpha);
end

Cm_ba = zeros(1, N); %coeficiente de momento de sustentacion de cada panel
for i=1: N 
    Cm_ba(i) = -Cl(i)*d(i)*cos(delta(i) - alpha);
end
Cm_ba_total = 0; %coeficiente de momento de sustentacion total
for i=1: N 
    Cm_ba_total = Cm_ba_total + long_panel(i)*Cm_ba(i);
end

Cm_ba1 = Cm_ba_total
x_cp = -Cm_ba_total/CL_total %Creo que es la posición de centro de presiones

%Primero sacar Cm_centro aerodinamico con la ultima formula, usando como punto generico A, el borde de ataque. Sacar x_ca, y
%luego se puede sacar el momento en cualquier punto

%para sacar el alpha que de Cl= 0, ir probando valores, y luego se podria
%usar regresion





























