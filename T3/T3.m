clear
clc
%Condições iniciais de operação

% Alimentacao (entrada)
F0=40;      % m3/h
E_oil_0=0.5;    % adimensional
E_s_0=0.001;    % adimensional
E_w_0=1-(E_oil_0 + E_s_0);      % adimensional
F_s_0=F0*E_s_0;

% Saida 
F_w=30;                    % m3/h
F_oil=10;

% Acúmulo
F_s=F_s_0;

% Massas específicas
p_oil=750;     % kg/m3
p_w=1000;      % kg/m3
p_s=1900;      % kg/m3

% Dimensões do separador
h_ant=8  ;    %m
A_sc=20  ; %m2

%altura dos componentes
h_oil=2; %m
h_w=4.99; %m
h_s=0.001; %m
h_ws=h_w+h_s;

%todo Alturas de SetPoint
h_ws_set=0.8*h_ant ;    %m

% Vazoes para quando o erro (variavel no setponi - variavel) for zero
F_ws_errozero = (E_w_0 + E_s_0)*F0;    % ft3/h

% Passo de tempo das equacoes diferenciais
dt=0.005;   % h

% Tempo de simulacao
%tf=20;       % h
%Altura máxima de areia
h_s_f=0.02;

% Valores iniciais dos controladores
intPID=0;
erro1V=0;
erro2V=erro1V;

intPIDFJ=0;
erro1T=0;
erro2T=erro1T;

% Constante do controlador do fluido refrigerante
Kj=-2;
Fdfj=0.1;
Fifj=0.2;

t_array = [];
h_w_array = [];
h_s_array = [];
h_oil_array = [];
h_ant_array = [];
F0_array = [];
F_oil_array = [];
F_w_array = [];

t=0;
 while (h_s<=h_s_f)

    % Imprimindo os resultados    
    if abs( eval(sym(t)/pi) - ceil(eval(sym(t)/pi)) ) <= 1e-3
       %fprintf(' Tempo=%0.3f   Ca=%0.4f   Cb=%0.4f   Cc=%0.4f   Cd=%0.4f   T=%0.3f   V=%0.3f   F=%0.3f   Tj=%0.3f   Fj=%0.3f\n', t,Ca,Cb,Cc,Cd,T,V,F,Tj,Fj );
    end
        
    % Controlador PID na Vazão
    %----------------------------------------------------------
    % A integral do PID na altura    
      intPIDFJ = intPIDFJ + (h_ws_set-h_ws)*dt;
    % A derivada do PID na altura       
      derPIDFJ = (erro2T-erro1T)/dt;
      erro1T = erro2T;
    % A parte proporcional do PID na altura
      proPIDFJ = (h_ws_set-h_ws);
    % A equação do PID na altura
      F_w = F_ws_errozero + Kj*proPIDFJ + Kj*Fdfj*derPIDFJ + (Kj/Fifj)*intPIDFJ;
    %----------------------------------------------------------   
 
    E_oil=E_oil_0 + 0.005*t;
    E_s=E_s_0;
    E_w=1-(E_oil+E_s);
    F_w=(F0*E_w)+F_s;
     
    % Calculando as derivadas
    dh_oil=(1/A_sc)*((F0*E_oil)-F_oil)*dt;
    dh_w=(1/A_sc)*((F0*E_w)-F_w)*dt;
    dh_s=(1/A_sc)*(F0*E_s)*dt;
    dh_ws=dh_w+dh_s;
    
    % Calculando as variáveis
    h_oil=h_oil+dh_oil;
    h_w=h_w+dh_w;
    h_s=h_s+dh_s;
    h_ws=h_w+h_s;
    
    %todo colocar o condicional de altura
    if ((h_oil+h_w+h_s)<h_ant)
        F_oil=0;
    else 
        F_oil=(F0*E_oil)+(dh_ws*A_sc);
    end 
    
    t_array = [t_array t];
    h_w_array = [h_w_array h_w+h_s];
    h_s_array = [h_s_array h_s];
    h_oil_array = [h_oil_array h_oil+h_w+h_s];
    h_ant_array = [h_ant_array h_ant];
    F0_array = [F0_array F0];
    F_oil_array = [F_oil_array F_oil];
    F_w_array = [F_w_array F_w];
    
    t=t+dt;
 end
 
    %Plotando os gráficos
    subplot(2,2,1); plot(t_array,F0_array);xlabel('tempo');ylabel('F0');                %hold on;    
    subplot(2,2,2); plot(t_array,F_oil_array); xlabel('tempo');ylabel('F_oil');       %hold on;
    subplot(2,2,3);
    plot(t_array,h_w_array,'b',t_array,h_s_array,'y',t_array, h_oil_array,'k',t_array,h_ant_array,'r'); xlabel('tempo');ylabel('Altura');            %hold on;
    subplot(2,2,4); plot(t_array,F_w_array); xlabel('tempo');ylabel('Vazão de água');             %hold on;     