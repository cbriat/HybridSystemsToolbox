clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A = [1 3;
     -1 2];

J = [0.5 0;
    0 0.5];

Ec = 0.75*eye(2);
Ed = 0.2*eye(2);

% Degree of the polynonmial variables in the SOS program
d = 6;

% Tentative maximum dwell-time
T(1) = 1e-3;
T(2) = 0.347;
echo off;

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[R,P,info] = RangeDT_Impulsive_MSStability(A,Ec,J,Ed,T,d);
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

z0 = vec(P);
Az = kron(A',eye(2))+kron(eye(2),A')+kron(Ec',Ec');
for(i=1:length(Tt))
    zT = expm(Az*Tt(i))*z0;
    XiT{i} = mat(zT);
    LMI{i} = J'*XiT{i}*J+Ed'*XiT{i}*Ed-P;
    vp(:,i) = eig(LMI{i});
end

disp('***************************************************')
disp('We find that The discrete-time LMI J''*Xi(theta)*J+Ed''*Xi(theta)*Ed-P <= lambda*I where lambda is')
disp(max(max(vp)))
disp('***************************************************')

%%
disp('Press any key to simulate the system');
pause 
seqT = (T(2)-T(1))*rand(20,1)+T(1);
options.domain = T(2);
E0 = [-5;5];
temp = 2*rand(2,2);
V0 = temp'*temp;

[t,E,V,h] = SimStochImpSyst(A,Ec,J,Ed,E0,V0,seqT,options);

x0 = [-5;5];
dt = 1e-3;
[tp,x,hp] = SimStochImpSystPaths(A,Ec,J,Ed,x0,dt,seqT,options);

%%%%%%

