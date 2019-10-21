clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A{1} = [-2 1;
        5 -3];

A{2} = [0.1 0;
        0.1 0.2];
N = size(A,2);

% Degree of the polynonmial variables in the SOS program
d = 4;

% Tentative range dwell-time
Tmin(1) = 1;
Tmin(2) = 0.001;

Tmax(1) = Inf;
Tmax(2) = 1.28;
echo off


% Press any key to solve for the conditions
disp('Press any key to solve for the SOS conditions.')
pause
% We solve the conditions
[R,P,info] = RangeDT_Switched_Stability(A,Tmin,Tmax,d);
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
mode = 1;
seqT = [];
Ncycle = 40;
Tmni = Tmax;
Tmni(isinf(Tmax))=0;
for(i=1:N)    
   Tm(i) = min([2*max(Tmni) Tmax(i)]); 
end
for(i=1:Ncycle)    
    seqT = [seqT;Tmin(mode)+rand*(Tm(mode)-Tmin(mode)) mode];
    mode = 3-mode;
end

x0 = [-5; 5];
options.domain(1) = Tmax(1);
options.domain(2) = Tmax(2);
[t,x,h] = SimSwSyst(A,x0,seqT,options);