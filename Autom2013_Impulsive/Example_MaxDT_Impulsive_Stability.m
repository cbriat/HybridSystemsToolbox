clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 3;
     -1 2];

J = 0.5*eye(2);

% Degree of the polynonmial variables in the SOS program
d = 6;

% Tentative maximum dwell-time
T = 0.4620;

% Press any key to solve for the conditions
echo off

disp('Press any key to solve for the SOS conditions.')
pause

% We solve the conditions
[R,P,info] = MaxDT_Impulsive_Stability(A,J,T,d);

disp('The problem has been solved sucessfully if info.pinf==0, info.dinf==0 and info.numerr==0.');
if(info.problem~=0)
   error('The problem has not been properly solved.'); 
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
pause

% We now double-check and verify that the associated hybrid conditions are satisfied
disp('***************************************************')
disp('The continuous-time LMI A''*P+P*A gives')
disp(A'*P+P*A)
disp('and its eigenvalues (must be positive) are')
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
seqT = T*rand(1,20);
x0 = [-5 5];
options = [];
[t,x,hf,hl] = SimImpSyst(A,J,x0,seqT,options);

