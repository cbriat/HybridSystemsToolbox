clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 0;
     1 2];

alpha = 0;
A = A + alpha*eye(size(A,1));
 
J = [1 1;
    1 3];

Bc = [1;
      0];

Bd = [0;
    0];

% Degree of the polynonmial variables in the SOS program
d = 1;

% Tentative minimum dwell-time
T = 0.1;


echo off
disp('Press any key to solve for the SOS conditions.')
pause

% We solve the conditions
[Kc,Kd,R,P,info] = MinDT_Impulsive_Stabilization(A,Bc,J,Bd,T,d);


disp('The problem has been solved sucessfully if info.pinf==0, info.dinf==0 and info.numerr==0.');
if(info.problem~=0)
   error('The problem has not been properly solved.'); 
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
pause

P = P/norm(P);
KcT = double(subs(Kc,'tau',T));
h = 1e-4;
Tt = [0:h:T];
Af = inline(A+Bc*Kc,'tau');

Phi = eye(size(A,1));
for(i=1:length(Tt))
    Phi = Phi + h*Af(Tt(i))*Phi;
    vp(:,i) = real(eig(Af(Tt(i))));
end

% We now double-check and verify that the associated hybrid conditions are satisfied
disp('***************************************************')
disp('The continuous-time LMI (A+Bc*Kc(T))''*P+P*(A+Bc*Kc(T)) gives')
disp((A+Bc*KcT).'*P+P*(A+Bc*KcT))
disp('and its eigenvalues (must be negative) are')
disp(eig((A+Bc*KcT)'*P+P*(A+Bc*KcT)))
disp('***************************************************')
disp(' ')
disp('***************************************************')
disp('The discrete-time LMI (J+Bd*Kd)''*Phi(T)''*P*Phi(T)*(J+Bd*Kd)-P gives')
disp((J+Bd*Kd)''*Phi''*P*Phi*(J+Bd*Kd)-P)
disp('and its eigenvalues (must be negative) are')
disp(eig((J+Bd*Kd)''*Phi''*P*Phi*(J+Bd*Kd)-P))
disp('***************************************************')

disp('Press any key to simulate the system');
pause
Kcf = inline(Kc,'tau');

A = A - alpha*eye(size(A,1));
seqT = T+3*T*rand(1,20);
x0 = [-1 1];
options.minDT = T;

% We simulate the system with the continuous-time output yc=Kc*x and the
% discrete-time output yd=Kd*x
[t,td,x,yc,yd,hf,hl] = SimImpSystOut(A+Bc*Kc,J+Bd*Kd,[],[1 1],x0,seqT,options);