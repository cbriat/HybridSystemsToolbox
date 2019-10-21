clear all 
close all
clc
echo on;

%% Matrices of the linear impulsive system
A{1} = [0 1;
        -5 1];

A{2} = [0 1;
        -1 5];
    
B{1} = 0*[0;1];
B{2} = [0;1];

N = 2;
    
% Degree of the polynonmial variables in the SOS program
d = 4;

% Tentative range dwell-time
Tmin(1) = 0.1;
Tmin(2) = 1;

Tmax(1) = 0.2;
Tmax(2) = 1.1;

% Press any key to solve for the conditions
echo off
pause
echo on
% We solve the conditions
[K,R,P,info] = RangeDT_Switched_Stabilization(A,B,Tmin,Tmax,d);
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
   Acl{i} = A{i}+B{i}*K{i};
end
for(i=1:Ncycle)    
    seqT = [seqT;Tmin(mode)+rand*(Tm(mode)-Tmin(mode)) mode];
    mode = 3-mode;
end

x0 = [-5; 5];
options.domain(1) = Tmax(1);
options.domain(2) = Tmax(2);
options.K = K;
[t,x,h,u] = SimSwSyst(Acl,x0,seqT,options);

% 
% h = 1e-4;
% for(i=1:N)
%    Tt{i} = [0:h:Tmax(i)];
%    Af = inline(A{i}+B{i}*K{i},'tau');
%    Kf = inline(K{i},'tau');
%    Kt{i} = [];
%    Phi{i} = eye(size(A{i},1));
%    for(k=1:length(Tt{i}))
%        Phi{i} = Phi{i} + h*Af(Tt{i}(k))*Phi{i};
%        Kt{i} = [Kt{i}; Kf(Tt{i}(k))];
%        vp(:,k) = real(eig(Af(Tt{i}(k))));
%    end
% end
