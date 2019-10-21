clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 1;
     1 -2];
alpha = 3;
A = A+alpha*eye(size(A,1));
 
Bc1 = [4;0];
Bc2 = 0*[1;0];

J = [3 1; 1 2];

Ec = [1 0;
    1 2];
Ed = 0.2*[1 0;
    1 -1];

Bd1 = [1;0];
Bd2 = 0*[0;0.1];

% Degree of the polynonmial variables in the SOS program
d = 1;

% Tentative maximum dwell-time
T = 0.1;
echo off

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[Kc,Kd,S,P,info] = MinDT_Impulsive_MSStabilization(A,Ec,Bc1,Bc2,J,Ed,Bd1,Bd2,T,d);
P = P/norm(P);

if(info.problem~=0)
   error('The problem has not been properly solved.'); 
else
   disp('The problem has been sucessfully solved since info.pinf==0, info.dinf==0 and info.numerr==0.');
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
P
pause

%% We check the conditions
A = A-alpha*eye(size(A,1));
h = 1e-4;
Tt = [0:h:T];
Af = inline(A+Bc1*Kc,'tau');
Ef = inline(Ec+Bc2*Kc,'tau');

X = P;
for(i=1:length(Tt))
    X = X + h*(Af(Tt(i))'*X+X*Af(Tt(i))+Ef(Tt(i))'*X*Ef(Tt(i)));
%     vp(:,i) = real(eig(Af(Tt(i))));
end

KcT = double(subs(Kc,'tau',T));

disp('***************************************************')
disp('The continuous-time LMI (A+Bc1*KcT)''*P+P*(A+Bc1*KcT)+(Ec+Bc2*KcT)''*P*(Ec+Bc2*KcT) gives')
disp((A+Bc1*KcT)''*P+P*(A+Bc1*KcT)+(Ec+Bc2*KcT)''*P*(Ec+Bc2*KcT))
disp('and its eigenvalues (must be positive) are')
disp(eig((A+Bc1*KcT)''*P+P*(A+Bc1*KcT)+(Ec+Bc2*KcT)''*P*(Ec+Bc2*KcT)))
disp('***************************************************')
disp(' ')
disp('***************************************************')
disp('The discrete-time LMI (J+Bd1*Kd)''*X*(J+Bd1*Kd)+(Ed+Bd2*Kd)''*X*(Ed+Bd2*Kd)-P gives')
disp((J+Bd1*Kd)''*X*(J+Bd1*Kd)+(Ed+Bd2*Kd)''*X*(Ed+Bd2*Kd)-P)
disp('and its eigenvalues (must be negative) are')
disp(eig((J+Bd1*Kd)''*X*(J+Bd1*Kd)+(Ed+Bd2*Kd)''*X*(Ed+Bd2*Kd)-P))
disp('***************************************************')

%%
disp('Press any key to simulate the system');
pause

seqT = T+T*rand(30,1)/2;
options.domain = T;
options.Kc = Kc;
options.Kd = Kd;
E0 = [-5;5];
temp = 2*rand(2,2)-1;
V0 = temp'*temp;

[t,E,V,h,uc,ud] = SimStochImpSyst(A+Bc1*Kc,Ec,J+Bd1*Kd,Ed,E0,V0,seqT,options);

x0 = [-5;5];
dt = 1e-3;
[tp,x,hp,ucp,udp] = SimStochImpSystPaths(A+Bc1*Kc,Ec,J,Ed+Bd1*Kd,x0,dt,seqT,options);
