clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [0 1;
     0 -1];

Ad = [0;1]*[-0.3410 -0.1332];
E1 = [0 0;0 0.1];

alpha = 0.1;
E2 = alpha*[0;1]*[-0.3410 -0.1332];

% Degree of the polynonmial variables in the SOS program
d = 2;

% Tentative maximum dwell-time
T(1) = 1e-3;
T(2) = 1;
echo off

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[S,P,info] = RangeDT_Sampled_Data_MSStability(A,Ad,E1,E2,T,d);
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
st = (T(2)-T(1))/201;
Tt = [T(1):st:T(2)];

n = size(A,1);
na = 2*n;
Abar = [A Ad;
        zeros(n,na)];

E1bar = [E1 zeros(n,n);
         zeros(n,na)];
         
E2bar = [zeros(n,n) E2;
         zeros(n,na)];
     
Jbar = [eye(n) zeros(n,n);
        eye(n) zeros(n,n)];

z0 = vec(P);
Az = kron(Abar',eye(na))+kron(eye(na),Abar')+kron(E1bar',E1bar')+kron(E2bar',E2bar');
for(i=1:length(Tt))
    zT = expm(Az*Tt(i))*z0;
    XiT{i} = mat(zT);
    LMI{i} = Jbar'*XiT{i}*Jbar-P;
    vp(:,i) = eig(LMI{i});
end
maxvp = max(max(vp));

disp('***************************************************')
disp('The maximum eigenvalue (must be negative) of the discrete-time LMI Jbar''*Phi(theta)''*P*Phi(theta)*Jbar-P for all theta is');
disp(maxvp)
disp('***************************************************')

disp('Press any key to simulate the system');
pause 

seqT = (T(2)-T(1))*rand(1,10)+T(1);
x0 = [-5 5 -5 5];
options = [];
%% FINISH THAT, THE ISSUE IS THAT THERE ARE TWO SOURCES OF NOISE
[t,td,x,yc,yd,hf,hl] = SimImpSystOut(Abar,J+Bd*Kd,0*Kc,Kd,x0,seqT,options);
