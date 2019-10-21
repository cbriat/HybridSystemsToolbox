clear all 
close all
clc
echo on;

%% Matrices of the linear switched system
A{1} = [0 1;
    -2 -1];

A{2} = [0 1;
    -9 -1];

N = 2;
% Degree of the polynonmial variables in the SOS program
d = 6;

% Tentative minimum dwell-time
T = 0.6222;
echo off

% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[R,P,info] = MinDT_Switched_Stability(A,T,d);
for(i=1:N)
    P{i} = P{i}/norm(P{i});
end

if(info.problem~=0)
   error('The problem has not been properly solved.'); 
else
   disp('The problem has been sucessfully solved since info.pinf==0, info.dinf==0 and info.numerr==0.');
end
disp('Press any key to verify the continuous- and discrete-time conditions with the obtained Lyapunov matrix');
for(i=1:N)
    txt = ['P' num2str(i) '='];
    disp(txt)
    disp(P{i})
end
pause

% We now double-check and verify that the associated hybrid conditions are satisfied
for(i=1:N)
   txt  = ['A' num2str(i) '''*P' num2str(i) '+P' num2str(i) '*A' num2str(i) ' ='];
   disp(txt)
   disp(A{i}'*P{i}+P{i}*A{i})
   disp('and has the eigenvalues')
   disp(eig(A{i}'*P{i}+P{i}*A{i}))    
end

for(i=1:N)
    for(j=1:N)
        if(i~=j)
            txt  = ['e^{A' num2str(i) '''*T}*P' num2str(i) '*e^{A' num2str(i) '*T}-P' num2str(j) '='];
            disp(txt)
            disp(expm(A{i}'*T)*P{i}*expm(A{i}*T)-P{j})
            disp('and has the eigenvalues')
            disp(eig(expm(A{i}'*T)*P{i}*expm(A{i}*T)-P{j}))
        end
    end   
end

disp('Press any key to simulate the system');
pause 
seqT = [T+rand(20,1) randi(2,20,1)];
x0 = [-5; 5];
options.domain(1) = T;
options.domain(2) = T;
[t,x,h] = SimSwSyst(A,x0,seqT,options);
