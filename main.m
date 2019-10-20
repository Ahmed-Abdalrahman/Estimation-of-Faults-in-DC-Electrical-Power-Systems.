%% Estimation of Faults in DC Electrical Power Systems.
% 
%   
% |.|
%  Student names: *Ahmed Abdalrahman Ahmed* and *Farid Farmani*
% |.|
% 
%
%  Course: ECE 602- Winter 2017

%% 
clc
clear

% The incidence matrix G 
G = zeros (18,18);

G(1,1)= 1;  G(2,1)=-1;   G(2,2)=1;  G(3,2)=-1; G(3,3)=1;
G(4,3)=-1;  G(4,4)=-1;   G(12,4)=1; G(4,5)=1;  G(5,5)=-1;
G(5,6)=1;   G(6,6)=-1;   G(6,7)=1;  G(7,7)=-1;
G(5,8)=1;   G(8,8)=-1;   G(8,9)=1;  G(9,9)=-1;

G(10,10)=1;  G(11,10)=-1; G(11,11)=1;  G(12,11)=-1; G(12,12)=1;
G(13,12)=-1; G(3,13)=1;   G(13,13)=-1; G(13,14)=1;  G(14,14)=-1;
G(14,15)=1;  G(15,15)=-1; G(15,16)=1;  G(16,16)=-1; G(14,17)=1;
G(17,17)=-1; G(17,18)=1;  G(18,18)=-1;

% The selection matrix SB that un-selects the boundary nodes of the circuit
% where KCL does not hold.
SB= ones (18,18); 

SB(1,1)=0; SB(10,10)=0; SB(9,9)=0; SB(7,7)=0; SB(16,16)=0; SB(18,18)=0;

% Ki (Diagonal Matrix) diagonal entry of KI could be a branch resistance 
%branch resistans
Res= [0.1 0 10^15 0 0 0 6 0 10 0.1 0 0 10^15 0 0 4 0 20];  
Ki = diag (Res);

% Kv unity Diagonal Matrix
Kv = eye(18);

% The Branches current vector J
J=  sym('J', [18 1]);

% Voltage drop on Branches vector V
V= sym('V', [18 1]);

% Nodal volatage vector E
E= sym('e', [18 1]);          

%The KCL equations for the currents j 
kcl= zeros (18,1);
% disp( sprintf( 'KCL Equations'));
SB.*G*J == kcl;

%The KVL relates voltage drops V  to the node voltages E 
% disp( sprintf( 'KVL Equations '));
G.'*E == V;

%Branch Constitutive Equations (BCE)
Fg = sym('Fg', [18 1]); 
% disp( sprintf( 'Branch Constitutive Equations (BCE)'));
- Ki*J + Kv*V == Fg;

%  
%%
% *-Observation equations*

Jms = sym('Jms', [18 1]);  % current sensors
Ems = sym('Ems', [18 1]); % Volatage sensors

Esrc = zeros (18,1); %  Voltage sources
Esrc(1,1)= 25.48;  Esrc(10,1)=24.83;

Egrd = zeros (18,1); %  ground nodes 
Egrd(7,1)= 0;  Egrd(9,1)= 0;  Egrd(16,1)=0;  Egrd(18,1)=0;

Fs = sym('Fs', [18 1]); % Fault offsets in voltage source
FV = sym('FV', [18 1]); % Fault offsets in voltage sensor
FI = sym('FI', [18 1]); % Fault offsets in Current sensor

SI = zeros (18,18); %  selection matrices for branches with current sensor
SI(1,1)= 1;  SI(5,5)= 1; SI(8,8)= 1;
SI(10,10)=1;  SI(14,14)=1; SI(17,17)=1;

SV = zeros (18,18); %  selection matrices for nodes with voltage sensor
SV(2,2)=1;   SV(3,3)=1;  SV(4,4)=1;  SV(5,5)=1; SV(8,8)=1;
SV(11,11)=1;  SV(12,12)=1; SV(13,13)=1; SV(14,14)=1; SV(17,17)=1;

Ss = zeros (18,18); %  selection matrices voltage sources
Ss(1,1)= 1;  Ss(10,10)= 1; 

Sg = zeros (18,18); %  selection matrices ground nodes
Sg(7,7)= 1;  Sg(9,9)= 1;  Sg(16,16)=1;  Sg(18,18)=1;

SV*Ems == SV*E + SV*FV;
SI*Jms == SI*J + FI;

Esrc == Ss*E + Fs;
Egrd == Sg*E;

%%
% *-STA model*

A = [SB.*G,       zeros(18,18), zeros(18,18); ...
    zeros(18,18),  eye(18),        -G.'     ; ... 
    Ki,              Kv,         zeros(18,18)];

X = [J;V;E];

F= [Fg; FI; FV; Fs];

B= [zeros(36,72); eye(18), zeros(18,54)];

% Y = [Jms;Ems;Esrc;Egrd];      % when faults are considered
Y = [zeros(36,1);Esrc;Egrd];    % no faults assumption


% C = [SI, zeros(18,36); ...    % when faults are considered
%     zeros(18,36), SV; ...
%     zeros(18,36), Ss; ...
%     zeros(18,36), Sg];

C = [ zeros(36,54); ...         % no faults assumption
    zeros(18,36), Ss; ...
    zeros(18,36), Sg];

% D = [zeros(18,18), SI, zeros(18,36); ...  % when faults are considered
%     zeros(18,36), SV, zeros(18,18); ...
%     zeros(18,54), Ss; zeros(18,72)];

D = [ zeros(36,72); ...        % no faults assumption
    zeros(18,54), Ss; zeros(18,72)];

% disp( sprintf( 'Circuit analysis equations '));
% A*X + B*F == 0              % Model equations
% Y == C*X + D*F




%% Visualization of results 
% The optimization problem (3) is solved to obtain the estimation of
% faults. two case studies are considered for this reason. These include
% base case with no faults in the system and two other cases with
% multiple faults in the system. 

%%
% |.|
%  _*Case (1)*: base case, no faults occurred_
%%
% |.|
%  In this case we only calculate the true states of the circuit based on
%  its parameters. Indeed, wrong data from sensors is ignored in this case.
%  
% 
%  
%  

cvx_begin 
variables XX(54,1) FF(72,1)
minimize  0.5 * (square_pos(norm((A*XX + B*FF),2)) + square_pos(norm((C*XX+D*FF - Y),2)))+ norm(FF,1)
cvx_end

% Display Results
JJ = XX(1:18,1);
VV = XX(19:36,1);
EE = XX(37:54,1);
FFg = FF(1:18,1);
FFI = FF(19:36,1);
FFV = FF(37:54,1);
FFs = FF(55:72,1);
disp( sprintf( '\nResults: Without any Faults\n=========================\ncvx_optval: %6.4f\ncvx_status: %s\n',cvx_optval, cvx_status ) );
disp( ['Branches Currents (A)'] );
disp( [ '   J  = [ ', sprintf( '%7.4f ', JJ.'), ']'] );
disp( [ sprintf( '\n'),'Branches Voltages (V)' ] ); 
disp( [ '   V  = [ ', sprintf( '%7.4f ', VV.'), ']'] );
disp( [ sprintf( '\n'),'Node Voltages (V)' ] ); 
disp( [ '   E  = [ ', sprintf( '%7.4f ', EE.'), ']'] );
disp( [sprintf( '\n'),'Faults in the branches'] );
disp( [ '   Fg  = [ ', sprintf( '%7.4f ', FFg.'), ']'] );
disp( [sprintf( '\n'), 'Fault offsets in Current' ] ); 
disp( [ '   Fi  = [ ', sprintf( '%7.4f ', FFI.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Voltages' ] ); 
disp( [ '   Fv  = [ ', sprintf( '%7.4f ', FFV.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Sources' ] ); 
disp( [ '   Fs  = [ ', sprintf( '%7.4f ', FFs.'), ']'] );


%% 
% _*Case (2)*: Optimization Problem with multiple faults in CB ISH280 and CS IT280_ 
% 
%
% In this case, we assume ISH280 fails to operate properly (opens) and current sensor IT280
% continues to show that it is closed. From the results, it can be seen
% that the model accurately estimate the faulted component and the value of
% the fault. In this case, the fault was in current sensor IT280 by
% (-0.1749 A), which is estimated in the Fault offsets in Current vector correctly.

KiErr =Ki;
KiErr(17,17)= 10^15; % circuit breaker resistance become Inf
Y3 = [JJ;EE;Esrc;Egrd]; % Sensors reading still as normal

AErr = [SB.*G,       zeros(18,18), zeros(18,18); ...  % when faults in CB is considered
    zeros(18,18),  eye(18),        -G.'     ; ... 
    KiErr,              Kv,         zeros(18,18)];

C = [SI, zeros(18,36); ...   % when faults are considered
    zeros(18,36), SV; ...
    zeros(18,36), Ss; ...
    zeros(18,36), Sg];

D = [zeros(18,18), SI, zeros(18,36); ... % when faults are considered
    zeros(18,36), SV, zeros(18,18); ...
    zeros(18,54), Ss; zeros(18,72)];

cvx_begin 
variables XX3(54,1) FF3(72,1)
minimize  0.5 * (square_pos(norm((AErr*XX3 + B*FF3),2)) + square_pos(norm((C*XX3+D*FF3 - Y3),2)))+ norm(FF3,1)
cvx_end

% Display Results
JJ3 = XX3(1:18,1);
VV3 = XX3(19:36,1);
EE3 = XX3(37:54,1);
FFg3 = FF3(1:18,1);
FFI3 = FF3(19:36,1);
FFV3 = FF3(37:54,1);
FFs3 = FF3(55:72,1);
disp( sprintf( '\n\nResults: With Faults in CB ISH280 and CS IT280' ) );
disp( sprintf( 'ISH280 fails (opens) and current sensor IT280 continues to show that it is closed.') );
disp( ['=========================================================='] );
disp( ['Branches Currents (A)'] );
disp( [ '   J  = [ ', sprintf( '%7.4f ', JJ3.'), ']'] );
disp( [ sprintf( '\n'),'Branches Voltages (V)' ] ); 
disp( [ '   V  = [ ', sprintf( '%7.4f ', VV3.'), ']'] );
disp( [ sprintf( '\n'),'Node Voltages (V)' ] ); 
disp( [ '   E  = [ ', sprintf( '%7.4f ', EE3.'), ']'] );
disp( [sprintf( '\n'),'Faults in the branches'] );
disp( [ '   Fg  = [ ', sprintf( '%7.4f ', FFg3.'), ']'] );
disp( [sprintf( '\n'), 'Fault offsets in Current' ] ); 
disp( [ '   Fi  = [ ', sprintf( '%7.4f ', FFI3.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Voltages' ] ); 
disp( [ '   Fv  = [ ', sprintf( '%7.4f ', FFV3.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Sources' ] ); 
disp( [ '   Fs  = [ ', sprintf( '%7.4f ', FFs3.'), ']'] );

%%
% _*Case (3)*: Optimization Problem with multiple faults in CB ISH180 and CS IT180_ 
% 
%
% In this case, we assume ISH180 fails to operate properly (opens) and current sensor IT180
% continues to show that it is closed. From the results, it can be seen
% that the model accurately estimate the faulted component and the value of
% the fault. In this case, the fault was in current sensor IT180 by
% (-1.3499 A), which is estimated in the Fault offsets in Current vector correctly.

KiErr =Ki;
KiErr(8,8)= 10^15; % circuit breaker resistance become Inf
Y5 = [JJ;EE;Esrc;Egrd]; % Sensors reading still as normal

AErr = [SB.*G,       zeros(18,18), zeros(18,18); ...  % when faults in CB is considered
    zeros(18,18),  eye(18),        -G.'     ; ... 
    KiErr,              Kv,         zeros(18,18)];

cvx_begin 
variables XX5(54,1) FF5(72,1)
minimize  0.5 * (square_pos(norm((AErr*XX5 + B*FF5),2)) + square_pos(norm((C*XX5+D*FF5 - Y5),2)))+ norm(FF5,1)
cvx_end

% Display Results
JJ5 = XX5(1:18,1);
VV5 = XX5(19:36,1);
EE5 = XX5(37:54,1);
FFg5 = FF5(1:18,1);
FFI5 = FF5(19:36,1);
FFV5 = FF5(37:54,1);
FFs5 = FF5(55:72,1);
disp( sprintf( '\n\nResults: With Faults in CB ISH180 and CS IT180' ) );
disp( sprintf( 'ISH180 fails (opens) and current sensor IT180 continues to show that it is closed.') );
disp( ['=========================================================='] );
disp( ['Branches Currents (A)'] );
disp( [ '   J  = [ ', sprintf( '%7.4f ', JJ5.'), ']'] );
disp( [ sprintf( '\n'),'Branches Voltages (V)' ] ); 
disp( [ '   V  = [ ', sprintf( '%7.4f ', VV5.'), ']'] );
disp( [ sprintf( '\n'),'Node Voltages (V)' ] ); 
disp( [ '   E  = [ ', sprintf( '%7.4f ', EE5.'), ']'] );
disp( [sprintf( '\n'),'Faults in the branches'] );
disp( [ '   Fg  = [ ', sprintf( '%7.4f ', FFg5.'), ']'] );
disp( [sprintf( '\n'), 'Fault offsets in Current' ] ); 
disp( [ '   Fi  = [ ', sprintf( '%7.4f ', FFI5.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Voltages' ] ); 
disp( [ '   Fv  = [ ', sprintf( '%7.4f ', FFV5.'), ']'] );
disp( [ sprintf( '\n'),'Fault offsets in Sources' ] ); 
disp( [ '   Fs  = [ ', sprintf( '%7.4f ', FFs5.'), ']'] );  

 

















