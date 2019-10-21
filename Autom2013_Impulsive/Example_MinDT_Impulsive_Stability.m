clear all 
close all
clc

echo on
%% Matrices of the linear impulsive system
A = [-1 0;
     1 -2];

J = [2 1;
    1 3];

% Degree of the polynonmial variables in the SOS program
d = 6;

% Tentative minimum dwell-time
T = 1.1406;
echo off

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[R,P,info] = MinDT_Impulsive_Stability(A,J,T,d);
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
disp('***************************************************')
disp('The continuous-time LMI A''*P+P*A gives')
disp(A'*P+P*A)
disp('and its eigenvalues (must be negative) are')
disp(eig(A'*P+P*A))
disp('***************************************************')
disp(' ')
disp('***************************************************')
disp('The discrete-time LMI J''*expm(A''*T)*P*expm(A*T)*J-P gives')
disp(J'*expm(A'*T)*P*expm(A*T)*J-P)
disp('and its eigenvalues (must be negative) are')
disp(eig(J'*expm(A'*T)*P*expm(A*T)*J-P))

%%
disp('Press any key to simulate the system');
pause 
seqT = T+rand(1,20);
x0 = [-5 5];
options.domain = T;
[t,x,hf,hl] = SimImpSyst(A,J,x0,seqT,options);