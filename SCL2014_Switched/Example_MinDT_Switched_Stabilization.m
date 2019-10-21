clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A{1} = [0 1;
        -5 1];

A{2} = [0 1;
        -1 5];
    
B{1} = [0;1];
B{2} = [0;1];

N = 2;
% Degree of the polynonmial variables in the SOS program
d = 1;

% Tentative minimum dwell-time
T = 0.1;
echo off


% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[K,R,P,info] = MinDT_Switched_Stabilization(A,B,T,d);
for(i=1:N)
    P{i} = P{i}/norm(P{i});
end

if(info.problem~=0)
   error('The problem has not been properly solved.'); 
else
   disp('The problem has been sucessfully solved since info.pinf==0, info.dinf==0 and info.numerr==0.');
end

disp('Press any key to simulate the system');
pause

for(i=1:N)
    KT{i} = double(subs(K{i},'tau',T));
end
for(i=1:N)
   Acl{i} = A{i}+B{i}*KT{i};    
end

seqT = [T+rand(20,1) randi(2,20,1)];
x0 = [-5; 5];
options.domain(1) = T;
options.domain(2) = T;
options.K = K;
for(i=1:N)
    Acl{i} = A{i}+B{i}*KT{i};
end
[t,x,h,u] = SimSwSyst(Acl,x0,seqT,options);
