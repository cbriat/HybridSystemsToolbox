clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 1;
     0 2];
 

alpha = 1;    % We can shift the matrix to ensure a certain decay rate
A = A + alpha*eye(size(A,1));

Bc = [0;1];

J = [1.2 1;
    0 1];

Bd = [1;0];

% Degree of the polynonmial variables in the SOS program
d = 2;

% Tentative range dwell-time
T(1) = 0.1;
T(2) = 0.5;

echo off
disp('Press any key to solve for the SOS conditions.')
pause

% We solve the conditions
[Kc,Kd,R,P,info] = RangeDT_Impulsive_Stabilization(A,Bc,J,Bd,T,d);
P = P/norm(P);

if(info.problem~=0)
   error('The problem has not been properly solved.'); 
else
   disp('The problem has been sucessfully solved since info.pinf==0, info.dinf==0 and info.numerr==0.');
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
P
pause

h = 1e-4;
Tt = [0:h:T(2)];
Af = inline(A+Bc*Kc,'tau');

Phi{1} = eye(size(A,1));
vp(1) = max(eig(Phi{1}));
cpt = 1;
for(i=1:length(Tt))
    Phi{i+1} = Phi{i} + h*Af(Tt(i))*Phi{i};
    if((Tt(i)>=T(1))&&(Tt(i)<=T(2)))
        vp(cpt) = max(eig((J+Bd*Kd)'*Phi{i+1}'*P*Phi{i+1}*(J+Bd*Kd)-P));
        cpt = cpt+1;
    end
end
maxvp = max(vp);

% We now double-check and verify that the associated discrete stability condition is satisfied
disp('***************************************************')
disp('The maximum eigenvalue (must be negative) of the discrete-time LMI (J+Bd*Kd)''*Phi(theta)''*P*Phi(theta)*(J+Bd*Kd)-P for all theta is');
disp(maxvp)

%%
disp('Press any key to simulate the system');
pause 
A = A - alpha*eye(size(A,1));
seqT = (T(2)-T(1))*rand(1,10)+T(1);
x0 = [-5 5];
options = [];

[t,td,x,yc,yd,hf,hl] = SimImpSystOut(A+Bc*Kc,J+Bd*Kd,0*Kc,Kd,x0,seqT,options);

 