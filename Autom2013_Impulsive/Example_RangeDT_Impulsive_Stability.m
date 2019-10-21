clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [-1 0.1;
     0 1.2];

J = blkdiag(1.2,0.5);

% Degree of the polynonmial variables in the SOS program
d = 8;

% Tentative range dwell-time
T(1) = 0.1824;
T(2) = 0.5776;

echo off
disp('Press any key to solve for the SOS conditions.')
pause

% We solve the conditions
[R,P,info] = RangeDT_Impulsive_Stability(A,J,T,d);


disp('The problem has been solved sucessfully if info.pinf==0, info.dinf==0 and info.numerr==0.');
if(info.problem~=0)
   error('The problem has not been properly solved.'); 
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
pause

% We now double-check and verify that the associated discrete stability condition is satisfied
syms theta
disp('***************************************************')
disp('The discrete-time LMI J''*expm(A''*theta)*P*expm(A*theta)*J-P gives')
disp(vpa(J'*expm(A'*theta)*P*expm(A*theta)*J-P,3))
disp('and the maximum of its eigenvalues (must be negative) for all theta is')
maxvp = -Inf;
st = (T(2)-T(1))/201;
Tt = [T(1):st:T(2)];
for(i=1:length(Tt))
    vp(:,i) = eig(expm(A'*Tt(i))*J'*P*J*expm(A*Tt(i))-P);
    maxvp = max([maxvp vp(:,i).']);
end
disp(maxvp);

%%
disp('Press any key to simulate the system');
pause 
seqT = (T(2)-T(1))*rand(1,30)+T(1);
x0 = [-5 5];
options.domain = T(2);
[t,x,hf,hl] = SimImpSyst(A,J,x0,seqT,options);

 