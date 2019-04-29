%Condições iniciais de operação

% Alimentacao (entrada)
F0=40;      % ft3/h
T0=530;     % °R
Ca0=0.5;    % lb mol/ft3
Cb0=0.5;    % lb mol/ft3
Cc0=0;      % lb mol/ft3
Cd0=0;      % lb mol/ft3
Fj=10.028;  % ft3/h
Tj0=530;    % °R

% Saida
F=40;                    % ft3/h
T=590.032;               % °R
Ca=0.3942;               % lb mol/ft3
Cb=0.3942;               % lb mol/ft3
Cc=0.1058;               % lb mol/ft3
Cd=0.1058;               % lb mol/ft3
Tj=588.400 ;              % °R

% Dimensoes
V=48;       % ft3
vj=3.85;    % ft3A=250;      
A=150;      % ft2

% Coeficiente global de troca térmica
U=150;      % BTU/(h ft2 °R)

% Propriedades do Fluido reativo
cp=0.75;    % BTU/(lb °R)
ro=50;      % lb/ft3

% Propriedades do Fluido refrigerante (da jaqueta)
cpj=1;      % BTU/(lb °R)
roj=62.3;   % lb/ft3

% Condicoes de operacao (set points)
Tset=590.032;   % °R
Vset=48;        % ft3

% Vazoes para quando o erro (variavel no setponi - variavel) for zero
Fjerrozero = 10.028;    % ft3/h
Ferrozero = 40;         % ft3/h

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

% Disturbio
%  F0=45;     % ft3/h

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

    dVdt=F0-F;
    dVCadt=F0*Ca0-F*Ca-V*k*Ca*Cb;
    dVCbdt=F0*Cb0-F*Cb-V*k*Ca*Cb;
    dVCcdt=F0*Cc0-F*Cc+V*k*Ca*Cb;
    dVCddt=F0*Cd0-F*Cd+V*k*Ca*Cb;
    dVTdt= F0*T0-F*T-deltaH*V*k*Ca*Cb/(cp*ro)-Q/(cp*ro);
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
    
    % Plotando os gráficos
    subplot(3,3,1); plot(tempo,Ca);xlabel('tempo');ylabel('Ca');                hold on;    
    subplot(3,3,2); plot(tempo,T); xlabel('tempo');ylabel('Temperatura');       hold on;
    subplot(3,3,3); plot(tempo,V); xlabel('tempo');ylabel('Volume');            hold on;
    subplot(3,3,4); plot(tempo,F); xlabel('tempo');ylabel('Vazão');             hold on;
    subplot(3,3,5); plot(tempo,Tj);xlabel('tempo');ylabel('Temp. Refigerante'); hold on;
    subplot(3,3,6); plot(tempo,Fj);xlabel('tempo');ylabel('Vazão Refrigerente');hold on;
    subplot(3,3,7); plot(tempo,Cb);xlabel('tempo');ylabel('Cb');                hold on;
    subplot(3,3,8); plot(tempo,Cc);xlabel('tempo');ylabel('Cc');                hold on;
    subplot(3,3,9); plot(tempo,Cd);xlabel('tempo');ylabel('Cd');                hold on;
    
    tempo = tempo + dt;
 end
