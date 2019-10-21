clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 1;
     1 2]/2;

J = [0.2 1;
    0 0.3];

Ec = 0.1*eye(2);
Ed = 0.1*eye(2);

% Degree of the polynonmial variables in the SOS program
d = 6;

% Tentative maximum dwell-time
T = 0.4;
echo off;

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[R,P,info] = MaxDT_Impulsive_MSStability(A,Ec,J,Ed,T,d);
P = P/norm(P);

if(info.problem~=0)
   error('The problem has not been properly solved.'); 
else
   disp('The problem has been sucessfully solved since info.pinf==0, info.dinf==0 and info.numerr==0.');
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
P
pause

% We now double-check and verify that the associated hybrid conditions are satisfied
z0 = vec(P);
Az = kron(A',eye(2))+kron(eye(2),A')+kron(Ec',Ec');
zT = expm(Az*T)*z0;
XiT = mat(zT);


disp('***************************************************')
disp('The continuous-time LMI A''*P+P*A+Ec''*P*Ec gives')
disp(A'*P+P*A+Ec'*P*Ec)
disp('and its eigenvalues (must be positive) are')
disp(eig(A'*P+P*A+Ec'*P*Ec))
disp('***************************************************')
disp(' ')
disp('***************************************************')
disp('The discrete-time LMI J''*XiT*J+Ed''*XiT*Ed-P gives')
disp(J''*XiT*J+Ed''*XiT*Ed-P)
disp('and its eigenvalues (must be negative) are')
disp(eig(J''*XiT*J+Ed''*XiT*Ed-P))

%%
disp('Press any key to simulate the system');
pause 
seqT = 0.99*T*rand(20,1);
options.domain = T;
E0 = [-5;5];
temp = 2*rand(2,2);
V0 = temp'*temp;

[t,E,V,h] = SimStochImpSyst(A,Ec,J,Ed,E0,V0,seqT,options);

x0 = [-5;5];
dt = 1e-3;
[tp,x,hp] = SimStochImpSystPaths(A,Ec,J,Ed,x0,dt,seqT,options);
