clear
clc
%Condições iniciais de operação

% Alimentacao (entrada)
F0=40;      % m3/h
E_oil_0=0.5;    % adimensional
E_s_0=0.01;    % adimensional
E_w_0=1-(E_oil_0 + E_s_0);      % adimensional

% Saida 
F_w=30;                    % m3/h
F_oil=10;

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
h_s=0.01; %m
h_ws=h_w+h_s;

%todo Alturas de SetPoint
h_ws_set=0.8*h_ant ;    %m

% Vazoes para quando o erro (variavel no setponi - variavel) for zero
F_ws_errozero = (E_w_0 + E_s_0)*F0;    % ft3/h

% Passo de tempo das equacoes diferenciais
dt=0.005;   % h

% Tempo de simulacao
tf=10;       % h

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


t=0;
 while (t<=tf)

    % Imprimindo os resultados    
    if abs( eval(sym(tempo)/pi) - ceil(eval(sym(tempo)/pi)) ) <= 1e-3
       fprintf(' Tempo=%0.3f   Ca=%0.4f   Cb=%0.4f   Cc=%0.4f   Cd=%0.4f   T=%0.3f   V=%0.3f   F=%0.3f   Tj=%0.3f   Fj=%0.3f\n', tempo,Ca,Cb,Cc,Cd,T,V,F,Tj,Fj );
    end
        
    % Controlador PID na Temperatura
    %----------------------------------------------------------
    % A integral do PID no volume    
      intPIDFJ = intPIDFJ + (h_ws_set-h_ws)*dt;
    % A derivada do PID no volume       
      derPIDFJ = (erro2T-erro1T)/dt;
      erro1T = erro2T;
    % A parte proporcional do PID no volume
      proPIDFJ = (h_ws_set-h_ws);
    % A equação do PID no volume
      F_w = F_ws_errozero + Kj*proPIDFJ + Kj*Fdfj*derPIDFJ + (Kj/Fifj)*intPIDFJ;
    %----------------------------------------------------------   
    
    
    % Controlador PID no Volume
    %----------------------------------------------------------
    % A integral do PID no volume    
      intPID = intPID + (Vset-V)*dt;
    % A derivada do PID no volume       
      derPID = (erro2V-erro1V)/dt;
      erro1V = erro2V;
    % A parte proporcional do PID no volume
      proPID = (Vset-V);
    % A equação do PID no volume
      Fw = Ferrozero + Kf*proPID + Kf*Tdf*derPID + (Kf/Tif)*intPID;
    %----------------------------------------------------------   
     
    E_oil=E_oil_0 + 0.1*sin(3*t);
    E_s=E_s_0;
    E_w=1-(E_oil+E_s);
    
    F_w=F0*E_w;
     
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
    
    % Plotando os gráficos
    %subplot(2,2,1); plot(t,F0);xlabel('tempo');ylabel('F0');                %hold on;    
    %subplot(2,2,2); plot(t,F_oil); xlabel('tempo');ylabel('F_oil');       %hold on;
    %subplot(2,2,3); plot(t,h_w); xlabel('tempo');ylabel('Volume');            %hold on;
    %subplot(2,2,4); plot(t,F_w); xlabel('tempo');ylabel('Vazão');             %hold on;
    
    
    t=t+dt;
 end