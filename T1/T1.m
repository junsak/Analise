clear
clc
R = 0.2;
%Condições iniciais de operação

% Alimentacao (entrada)
F0=40;      % ft3/h
T0=670;     % °R
Ca0=0.5;    % lb mol/ft3
Cb0=0.5;    % lb mol/ft3
Cc0=0;      % lb mol/ft3
Cd0=0;      % lb mol/ft3
Fj=24.527;  % ft3/h
Tj0=530;    % °R

% Saida
F=50;                    % ft3/h
T=670;               % °R
Ca=0.155;               % lb mol/ft3
Cb=0.155;               % lb mol/ft3
Cc=0.307;               % lb mol/ft3
Cd=0.307;               % lb mol/ft3
Tj=658.959;              % °R
CcF=Cc*F;

% Refluxo de saida (Nó 2) 
F2=F*R;
T2=T;
Ca2=Ca;
Cb2=Cb;
Cc2=Cc;
Cd2=Cd;

% Nó 1
F1=F2+F0;
T1=((F0*T0)+(F2*T2))/(F0+F2);
Ca1=((F0*Ca0)+(F2*Ca2))/(F0+F2);
Cb1=((F0*Cb0)+(F2*Ca2))/(F0+F2);
Cc1=((F0*Cc0)+(F2*Ca2))/(F0+F2);
Cd1=((F0*Cd0)+(F2*Ca2))/(F0+F2);

% Dimensoes
V=48;       % ft3
vj=3.85;    % ft3A=250;      
A=150;      % ft2

% Coeficiente global de troca térmica
U=250;      % BTU/(h ft2 °R)

% Propriedades do Fluido reativo
cp=0.75;    % BTU/(lb °R)
ro=50;      % lb/ft3

% Propriedades do Fluido refrigerante (da jaqueta)
cpj=1;      % BTU/(lb °R)
roj=62.3;   % lb/ft3

% Condicoes de operacao (set points)
Tset=670;   % °R
Vset=48;        % ft3

% Vazoes para quando o erro (variavel no setponi - variavel) for zero
Fjerrozero = 51.534;    % ft3/h
Ferrozero = 50;         % ft3/h

% Calor de reacao
deltaH=-30000;          %BTU/lb mol de A

% Constantes de Arrenius
Aa=7.08e10; % fator pre-exponencial  1/h
Ea=30000;   % energia de ativação    BTU/lb mol
Ra=1.99;    % constante dos gases    BTU/(lb mol °R)

% Constante do controlador do fluido refrigerante
Kj=-2;
Tdfj=0.1;
Tifj=0.2;

% Constantes do controlador da vazao de saida
Kf=-5;
Tdf=0.1;
Tif=0.2;

% Passo de tempo das equacoes diferenciais
dt=0.005;   % h

% Passo de impressão da tabela
pi=0.2;

% Tempo de simulacao
tf=10;       % h

% Calculado os produtos V*Ca, V*Cb e V*T
VCa=V*Ca;
VCb=V*Cb;
VCc=V*Cc;
VCd=V*Cd;
VT=V*T;

% Valores iniciais dos controladores
intPID=0;
erro1V=0;
erro2V=erro1V;

intPIDFJ=0;
erro1T=0;
erro2T=erro1T;

% Definição das matrizes das variáveis

Ca_array = [];
Cb_array = [];
Cc_array = [];
Cd_array = [];
F_array = [];
T_array = [];
Fj_array = [];
Tj_array = [];
V_array = [];
tempo_array = [];
CcF_array = [];

% Disturbio
%F0=45;     % ft3/h

 % Inicio do loop
 tempo=0;
 while (tempo<=tf)

    % Imprimindo os resultados    
    if abs( eval(sym(tempo)/pi) - ceil(eval(sym(tempo)/pi)) ) <= 1e-3
       fprintf(' Tempo=%0.3f   Ca=%0.4f   Cb=%0.4f   Cc=%0.4f   Cd=%0.4f   T=%0.3f   V=%0.3f   F=%0.3f   Tj=%0.3f   Fj=%0.3f\n', tempo,Ca,Cb,Cc,Cd,T,V,F,Tj,Fj );
    end
        
    % Controlador PID na Temperatura
    %----------------------------------------------------------
    % A integral do PID no volume    
      intPIDFJ = intPIDFJ + (Tset-T)*dt;
    % A derivada do PID no volume       
      derPIDFJ = (erro2T-erro1T)/dt;
      erro1T = erro2T;
    % A parte proporcional do PID no volume
      proPIDFJ = (Tset-T);
    % A equação do PID no volume
      Fj = Fjerrozero + Kj*proPIDFJ + Kj*Tdfj*derPIDFJ + (Kj/Tifj)*intPIDFJ;
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
      F = Ferrozero + Kf*proPID + Kf*Tdf*derPID + (Kf/Tif)*intPID;
    %----------------------------------------------------------   

    % k da reação
    k=(Aa)*exp(-Ea/(Ra*T));

    % O calor trocado com a camisa
    Q=U*A*(T-Tj);

    % Calculando as derivadas

    dVdt=F1-F;
    dVCadt=F1*Ca1-F*Ca-V*k*Ca*Cb;
    dVCbdt=F1*Cb1-F*Cb-V*k*Ca*Cb;
    dVCcdt=F1*Cc1-F*Cc+V*k*Ca*Cb;
    dVCddt=F1*Cd1-F*Cd+V*k*Ca*Cb;
    dVTdt= F1*T1-F*T-deltaH*V*k*Ca*Cb/(cp*ro)-Q/(cp*ro);
    dTjdt=Fj*(Tj0-Tj)/vj + Q/(cpj*roj*vj);

    % Calculando as variáveis
    V=V+dVdt*dt;
    VCa=VCa+dVCadt*dt;
    VCb=VCb+dVCbdt*dt;
    VCc=VCc+dVCcdt*dt;
    VCd=VCd+dVCddt*dt;    
    VT=VT+dVTdt*dt;
    Tj=Tj+dTjdt*dt;
    Ca=VCa/V;
    Cb=VCb/V;
    Cc=VCc/V;
    Cd=VCd/V;    
    T=VT/V;
    
    erro2V = Vset-V;
    erro2T = Tset-T;
    
    tempo = tempo + dt;
    
        %Saída
    Ts=T;
    Fs=F-F2;
    Cas=Ca;
    Cbs=Cb;
    Ccs=Cc;
    Cds=Cd;
    
    %Recalculando o reciclo
    F2=F*R;
    T2=T;
    Ca2=Ca;
    Cb2=Cb;
    Cc2=Cc;
    Cd2=Cd;

    %Recalcular a entrada 1 com novo reciclo(Nó 1)
    F1=F0+F2;
    T1=((F0*T0)+(F2*T2))/(F0+F2);
    Ca1=((F0*Ca0)+(F2*Ca2))/(F0+F2);
    Cb1=((F0*Cb0)+(F2*Ca2))/(F0+F2);
    Cc1=((F0*Cc0)+(F2*Ca2))/(F0+F2);
    Cd1=((F0*Cd0)+(F2*Ca2))/(F0+F2);
    
    Ca_array = [Ca_array Ca];
    Cb_array = [Cb_array Cb];
    Cc_array = [Cc_array Cc];
    Cd_array = [Cd_array Cd];
    F_array = [F_array F];
    T_array = [T_array T];
    Fj_array = [Fj_array Fj];
    Tj_array = [Tj_array Tj];
    V_array = [V_array V];
    tempo_array = [tempo_array tempo];
    CcF_array = [CcF_array CcF];
    
 end

     % Plotando os gráficos
    subplot(3,4,1); plot(tempo_array,Ca_array,'.');xlabel('tempo');ylabel('Ca');             hold on;    
    subplot(3,4,2); plot(tempo_array,T_array, '.'); xlabel('tempo');ylabel('Temperatura');   hold on;
    subplot(3,4,3); plot(tempo_array,V_array, '.'); xlabel('tempo');ylabel('Volume');            hold on;
    subplot(3,4,4); plot(tempo_array,F_array, '.'); xlabel('tempo');ylabel('Vazão');             hold on;
    subplot(3,4,5); plot(tempo_array,Tj_array, '.');xlabel('tempo');ylabel('Temp. Refigerante'); hold on;
    subplot(3,4,6); plot(tempo_array,Fj_array, '.');xlabel('tempo');ylabel('Vazão Refrigerente');hold on;
    subplot(3,4,7); plot(tempo_array,Cb_array, '.');xlabel('tempo');ylabel('Cb');                hold on;
    subplot(3,4,8); plot(tempo_array,Cc_array, '.');xlabel('tempo');ylabel('Cc');                hold on;
    subplot(3,4,9); plot(tempo_array,Cd_array, '.');xlabel('tempo');ylabel('Cd');                hold on;
    subplot(3,4,10); plot(tempo_array,CcF_array, '.');xlabel('tempo');ylabel('Cc*F');            hold on;
    
    drawnow