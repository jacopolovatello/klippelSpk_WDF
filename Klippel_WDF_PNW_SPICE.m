%% Jacopo Lovatello, M.862780 CP. 10546576
%% Progetto SASP A.A. 2015/2016

clear, close all, clc; 

comp = 1;           % 1 = Vout, 2 = Vin, 3 = Vtrasf  

%% Signal Variables Declaration 
Fs = 48000;         % Sampling Freq. 
N = Fs/100;         % Samples for 0.01 seconds execution
t = 0:N-1;          % Time axis 

f = 1000;           % Input Freq.
bias = 0;           % DC Bias
gain = 1;           % Gain

u = bias - gain*sin((2*pi*f/Fs).*t); 

%% Circuit Variables declaration 
Re = 3.02;        
Le = 0.132e-3; 
L2m = 0.264e-3; 
R2m = 3.52; 
Cms = 0.43e-3;    % mm/N
Mms = 11.87e-3;   % g 
Rms = 2.27;       % kg/sec 
Bl = 5.23;        % N/A

% Wave variables init 
a = zeros(N, 15);           % Null C.I.   INs  
b = zeros(N, 15);           % Null C.I.   OUTs

Vout = []; 
Vout_PNW = []; 

% Adaptors' ports definition according to WAVE scheme 
R1 = Re;                    % Adapted linear resistance, real generator  
R2 = 2*Le*Fs;               % Adapted inductor
R3 = R1 + R2;               % Series adaptor 
R4 = R3;                    
R7 = R2m;                   % Adapted linear resistance 
R8 = 2*L2m*Fs;              % Adapted inductor 
R6 = (R7*R8)/(R7+R8);       % Parallel adaptor 
R5 = R6; 
R9 = R5 + R4;               % Series adaptor 
R10 = R9*(Bl^2);            % Transformer 
R11 = 1/(2*Cms*Fs);         % Adapted capacitor 
R12 = R10 + R11;            % Series adaptor 
R13 = R12; 
R14 = 2*Mms*Fs;             % Adapted inductor 
R15 = R13 + R14;            % Series adaptor 

%% Wave implementation 
% Series adaptors' coefficients 
gamma1 = R1/(R1+R2); 
gamma2 = R4/(R4+R5); 
gamma3 = R10/(R10+R11); 
gamma4 = R13/(R13+R14); 

series123 =     [(1-gamma1) -gamma1 -gamma1; (-1+gamma1) gamma1 (-1+gamma1); -1 -1 0]; 
series459 =     [(1-gamma2) -gamma2 -gamma2; (-1+gamma2) gamma2 (-1+gamma2); -1 -1 0];
series101112 =  [(1-gamma3) -gamma3 -gamma3; (-1+gamma3) gamma3 (-1+gamma3); -1 -1 0]; 
series131415 =  [(1-gamma4) -gamma4 -gamma4; (-1+gamma4) gamma4 (-1+gamma4); -1 -1 0];

% Parallel adaptors' coefficients 
alfa5 = (1/R7)/((1/R7)+(1/R8));
parallel786 = [(alfa5-1) (1-alfa5) 1; alfa5 -alfa5 1; alfa5 (1-alfa5) 0]; 

% Non-adapted linear resistance 
r15 = (Rms-R15)/(Rms+R15); 

% Transformer normalized PW coefficient 
Bl_PNW = 1;

%% Routine waves 
for ii = 2:N
    
    a(ii, 1) = u(ii);                                           % Real voltage source 
    a(ii, 2) = -b(ii-1, 2);                                     % Inductor 
    
    A = [a(ii, 1); a(ii, 2); a(ii, 3)]; 
    b(ii, 3) = series123(3,:)*A;                                % Series adaptor
    a(ii, 4) = b(ii, 3); 
    
    a(ii, 7) = 0;                                               % Linear resistor 
    
    a(ii, 8) = -b(ii-1, 8);                                     % Inductor
    
    A = [a(ii, 7); a(ii, 8); a(ii, 6)]; 
    b(ii, 6) = parallel786(3,:)*A;                              % Parallel adaptor 
    a(ii, 5) = b(ii, 6); 
    
    A = [a(ii, 4); a(ii, 5); 0]; 
    b(ii, 9) = series459(3,:)*A;                                % Series adaptor 
    
    a(ii, 10) = b(ii, 9)*(-Bl);                                 % Transformer 
    
    a(ii, 11) = b(ii-1, 11);                                    % Inductor 
    
    A = [a(ii, 10); a(ii, 11); 0]; 
    b(ii, 12) = series101112(3,:)*A;                            % Series adaptor 
    a(ii, 13) = b(ii, 12); 
    
    a(ii, 14) = -b(ii-1, 14);                                   % Inductor
    
    A = [a(ii, 13); a(ii, 14); a(ii, 15)]; 
    b(ii, 15) = series131415(3,:)*A;                            % Series adaptor
    
    a(ii, 15) = b(ii, 15) * r15;                                % Linear resistor NOT ADAPTED 
    
    A = [a(ii, 13); a(ii, 14); a(ii, 15)]; 
    b(ii, 13) = series131415(1,:)*A;                            % Series adaptor
    b(ii, 14) = series131415(2,:)*A;                            % Series adaptor  

    a(ii, 12) = b(ii, 13); 
    
    A = [a(ii, 10); a(ii, 11); a(ii, 12)]; 
    b(ii, 10) = series101112(1,:)*A;                            % Series adaptor
    b(ii, 11) = series101112(2,:)*A;                            % Series adaptor
    
    a(ii, 9) = b(ii, 10)/(-Bl);                                 % Transformer 
    
    A = [a(ii, 4); a(ii, 5); a(ii, 9)]; 
    b(ii, 4) = series459(1,:)*A;                                % Series adaptor
    b(ii, 5) = series459(2,:)*A;                                % Series adaptor
    
    a(ii, 3) = b(ii, 4); 
    a(ii, 6) = b(ii, 5); 
    
    A = [a(ii, 7); a(ii, 8); a(ii, 6)]; 
    b(ii, 7) = parallel786(1,:)*A;                              % Parallel adaptor
    b(ii, 8) = parallel786(2,:)*A;                              % Parallel adaptor
    
    A = [a(ii, 1); a(ii, 2); a(ii, 3)]; 
    b(ii, 1) = series123(1,:)*A;                                % Series adaptor
    b(ii, 2) = series123(2,:)*A;                                % Series adaptor
    
    if comp == 1
        Vout = [Vout ((a(ii, 15)+b(ii, 15))/2)];
    end 
    if comp == 2
        Vout = [Vout ((a(ii, 1)+b(ii, 1))/2)];
    end 
    if comp == 3
        Vout = [Vout ((a(ii, 10)+b(ii, 10))/2)];
    end 
    
end 

%% PNW implementation 
% Variables reset for PNW algorithm 
a = zeros(N, 15);           % Null C.I.   INs  
b = zeros(N, 15);           % Null C.I.   OUTs

% Series adaptors' coefficients 
gamma1 = (2*[R1,R2,R3]')/(R1+R2+R3);  
gamma2 = (2*[R4,R5,R9]')/(R4+R5+R9); 
gamma3 = (2*[R10,R11,R12]')/(R10+R11+R12);
gamma4 = (2*[R13,R14,R15]')/(R13+R14+R15);

series123 =     eye(3) - sqrt(gamma1)*sqrt(gamma1)'; 
series459 =     eye(3) - sqrt(gamma2)*sqrt(gamma2)'; 
series101112 =  eye(3) - sqrt(gamma3)*sqrt(gamma3)'; 
series131415 =  eye(3) - sqrt(gamma4)*sqrt(gamma4)'; 

% Parallel adaptors' coefficients 
alfa5 = (2*[1/R7,1/R8,1/R6]'/((1/R7)+(1/R8)+(1/R6))); 
parallel786 = -eye(3)+sqrt(alfa5)*sqrt(alfa5)'; 

%% Routine PNW 
for ii = 2:N
    
    a(ii, 1) = u(ii)/(2*sqrt(R1));                              % Real voltage source 
    a(ii, 2) = -b(ii-1, 2);                                     % Inductor 
    
    A = [a(ii, 1); a(ii, 2); a(ii, 3)]; 
    b(ii, 3) = series123(3,:)*A;                                % Series adaptor
    a(ii, 4) = b(ii, 3); 
    
    a(ii, 7) = 0;                                               % Linear resistor 
    
    a(ii, 8) = -b(ii-1, 8);                                     % Inductor
    
    A = [a(ii, 7); a(ii, 8); a(ii, 6)]; 
    b(ii, 6) = parallel786(3,:)*A;                              % Parallel adaptor 
    a(ii, 5) = b(ii, 6); 
    
    A = [a(ii, 4); a(ii, 5); 0]; 
    b(ii, 9) = series459(3,:)*A;                                % Series adaptor 
    
    a(ii, 10) = b(ii, 9)*(-Bl_PNW);                             % Transformer 
    
    a(ii, 11) = b(ii-1, 11);                                    % Inductor 
    
    A = [a(ii, 10); a(ii, 11); 0]; 
    b(ii, 12) = series101112(3,:)*A;                            % Series adaptor 
    a(ii, 13) = b(ii, 12); 
    
    a(ii, 14) = -b(ii-1, 14);                                   % Inductor
    
    A = [a(ii, 13); a(ii, 14); a(ii, 15)]; 
    b(ii, 15) = series131415(3,:)*A;                            % Series adaptor
    
    a(ii, 15) = b(ii, 15) * r15;                                % Linear resistor NOT ADAPTED 
    
    A = [a(ii, 13); a(ii, 14); a(ii, 15)]; 
    b(ii, 13) = series131415(1,:)*A;                            % Series adaptor
    b(ii, 14) = series131415(2,:)*A;                            % Series adaptor  

    a(ii, 12) = b(ii, 13); 
    
    A = [a(ii, 10); a(ii, 11); a(ii, 12)]; 
    b(ii, 10) = series101112(1,:)*A;                            % Series adaptor
    b(ii, 11) = series101112(2,:)*A;                            % Series adaptor
    
    a(ii, 9) = b(ii, 10)/(-Bl_PNW);                             % Transformer 
    
    A = [a(ii, 4); a(ii, 5); a(ii, 9)]; 
    b(ii, 4) = series459(1,:)*A;                                % Series adaptor
    b(ii, 5) = series459(2,:)*A;                                % Series adaptor
    
    a(ii, 3) = b(ii, 4); 
    a(ii, 6) = b(ii, 5); 
    
    A = [a(ii, 7); a(ii, 8); a(ii, 6)]; 
    b(ii, 7) = parallel786(1,:)*A;                              % Parallel adaptor
    b(ii, 8) = parallel786(2,:)*A;                              % Parallel adaptor
    
    A = [a(ii, 1); a(ii, 2); a(ii, 3)]; 
    b(ii, 1) = series123(1,:)*A;                                % Series adaptor
    b(ii, 2) = series123(2,:)*A;                                % Series adaptor
    
    if comp == 1
        Vout_PNW = [Vout_PNW (((a(ii, 15)+b(ii, 15))*sqrt(R15)))];
    end
    if comp == 2
        Vout_PNW = [Vout_PNW (((a(ii, 1)+b(ii, 1))*sqrt(R1)))];
    end
    if comp == 3
        Vout_PNW = [Vout_PNW (((a(ii, 10)+b(ii, 10))*sqrt(R10)))];
    end
        
end

%% SPICE data import
filename = 'Klippel_Circuit.txt'; 
data = tdfread(filename, '\t'); 

SPICE_time = data.time; 
if comp == 1
    SPICE_Vout = data.V0x28n0060x29;
end 
if comp == 2
    SPICE_Vout = data.V0x28n0030x29;
end 
if comp == 3
    SPICE_Vout = data.V0x28n0040x29;
end 


% Costant resample 
xcoord = (1:N)/Fs;
RESAMP_time = xcoord; 
RESAMP_Vout = interp1(SPICE_time, SPICE_Vout, RESAMP_time, 'linear');

%% Plot
figure(1)
subplot(2,1,1)
plot(1:N-1, Vout, 'b', 'LineWidth', 2), hold on 
plot(1:N-1, Vout_PNW, '--g', 'LineWidth', 2), hold on 
plot(1:N-1, RESAMP_Vout(1:N-1),'r', 'LineWidth', 2), hold off; 
legend('WAVE', 'PNW', 'SPICE'); 
axis tight;
title('WAVE-PNW-SPICE comparison, 1KHz sin'); 
xlabel('samples','interpreter','latex');
ylabel('$V_{out}$ (V)','interpreter','latex');

subplot(2,1,2)
plot(1:N-1, abs(RESAMP_Vout(1:N-1) - Vout), 'r', 'LineWidth', 2), hold on 
plot(1:N-1, abs(Vout_PNW - Vout), 'b', 'LineWidth', 2), hold on 
legend('SPICE-WAVE', 'PNW-WAVE'); 
axis tight ;
xlabel('samples','interpreter','latex');
ylabel('$\Delta$ (V)','interpreter','latex');

