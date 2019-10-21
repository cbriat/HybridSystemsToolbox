clear all 
close all
clc
echo on;

%% Matrices of the sampled-data system
A = [0 1;
     0 -0.1];
 
alpha = 1;
A = A+alpha*eye(size(A,1));

B = [0;0.1];
n = size(A,1);
m = size(B,2);

% Degree of the polynonmial variables in the SOS program
d = 1;

% Tentative range dwell-time
T(1) = 1e-2;
T(2) = 0.1;

echo off
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[K,S,P,info] = RangeDT_Sampled_Data_Stabilization(A,B,T,d);

disp('The problem has been solved sucessfully if info.pinf==0, info.dinf==0 and info.numerr==0.');
if(info.problem~=0)
   error('The problem has not been properly solved.'); 
end
disp('Press any key to verify the discrete-time condition with the obtained Lyapunov matrix');
pause


Aa = [A B;
      zeros(m,n+m)];
Ja = [eye(n) zeros(n,m);
      K];
syms theta
Phi = expm(Aa*theta)*Ja;
% We normalize the matrix to get a reasonable magnitude
P = P/norm(P);

maxvp = -Inf;
st = (T(2)-T(1))/201;
Tt = [T(1):st:T(2)];
M = inline(Phi);
for(i=1:length(Tt))
    vp(:,i) = eig(M(Tt(i))'*P*M(Tt(i))*-P);
    maxvp = max([maxvp vp(:,i).']);
end


% We now double-check and verify that the associated discrete stability condition is satisfied
disp('***************************************************')
disp('The discrete-time LMI Ja''*expm(Aa''*theta)*P*expm(Aa*theta)*Ja-Pa gives')
disp(vpa(Phi.'*P*Phi-P,3))
disp('and the maximum of its eigenvalues for all theta is')
disp(maxvp)
disp('***************************************************')


disp('Press any key to simulate the system');
pause 

A = A-alpha*eye(size(A,1));
Aa = [A B;
      zeros(m,n+m)];
seqT = (T(2)-T(1))*rand(1,150)+T(1);
x0 = [-1 1];
      
[t,x,hf,hl] = SimImpSyst(Aa,Ja,[x0 K*[x0 0]'],seqT,options);

% [t,td,x,yc,yd,hf,hl] = SimImpSystOut(Aa,Ja,zeros(1,n+m),K,[x0 K*[x0 0]'],seqT,options);