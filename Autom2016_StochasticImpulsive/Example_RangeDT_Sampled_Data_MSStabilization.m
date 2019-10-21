clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [0 1;
     0 -1];
 
alpha = 0;
A = A + alpha*eye(size(A,1));
 
B = [0;
    1];

E1 = [0 0;
      0 0.1];

alpha = 0.1;
E2 = alpha*B;

% Degree of the polynonmial variables in the SOS program
d = 1;

% Tentative maximum dwell-time
T(1) = 1e-3;
T(2) = 1;
echo off

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[K,S,P,info] = RangeDT_Sampled_Data_MSStabilization(A,B,E1,E2,T,d);
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
m = size(B,2);
na = n+m;
A = A - alpha*eye(size(A,1));

Abar = [A B;
        zeros(m,na)];

E1bar = [E1 zeros(n,m);
         zeros(m,na)];
         
E2bar = [zeros(n,n) E2;
         zeros(m,na)];
     
Jbar = [eye(n) zeros(n,m);
        K];

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


%%%% MULTIPLE SOUURCES OF NOISE CHECK THAT
% % seqT = (T(2)-T(1))*rand(30,1)+T(1);
% % options.domain = T(2);
% % E0 = [-5;5];
% % temp = 2*rand(2,2)-1;
% % V0 = temp'*temp;
% % 
% % [t,E,V,h] = SimStochImpSyst(A+Bc1*Kc,Ec,J+Bd1*Kd,Ed,E0,V0,seqT,options);
% % 
% % x0 = [-5;5];
% % dt = 1e-3;
% % [tp,x,hp] = SimStochImpSystPaths(Abar,,J,Ed+Bd1*Kd,x0,dt,seqT,options);
