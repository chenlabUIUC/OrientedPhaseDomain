% NAME:     CalculateDomainCenterofMass
% PURPOSE:  This code calculates the center of mass for each domain in a
%           cathode particle
% INPUT:
%           Domain map: 'DomainMap.xlsx'
% OUTPUT:
%           Center of mass for each domain: 'DomainCenterofMassCoordinates.xlsx'
%               1st column: x cordinates of [100]t domain
%               2nd column: y cordinates of [100]t domain
%               3rd column: x cordinates of [111]t domain
%               4th column: y cordinates of [111]t domain
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019

% load data
filename = 'DomainMap.xlsx';
B = xlsread(filename);
C=B;
D=B;
[XX,YY] = size(B);

% Select the [100]t and [111]t domain distributions, respectively
for x=1:1:XX
    for y=1:1:YY
        if B(x,y)~=2
            C(x,y)=0; % [100]t domain distribution
        end
        if B(x,y)~=1
            D(x,y)=0; % [111]t domain distribution
        end
        
    end
end

% label the each [100]t (or [111]t) domain with indices. 0 index is usually for
% background. Domain index starts from 1.
L = bwlabel(C,4);
M = bwlabel(D,4);

% Optional: show the labelled domains
%{
figure
P=imagesc(L);
hold on
figure
P=imagesc(M);
hold off
%}

% calculate the center of mass coordinates for [100]t domains in the particle
N1=max(max(L));
Result=[0,0];
for n=1:1:N1
    count=0;
    Corx=0;
    Cory=0;

    for x=1:1:XX
    for y=1:1:YY
        if L(x,y)==n
            Corx=Corx+x; 
            Cory=Cory+y;
            count=count+1;
        end
    end
    end
    
    Result(n,1)=Corx/count;  % x coordinate for the domain
    Result(n,2)=Cory/count;  % y coordinate for the domain
end

% calculate the center of mass coordinates for [111]t domains in the particle
N2=max(max(M));
for n=1:1:N2
    count=0;
    Corx=0;
    Cory=0;

    for x=1:1:XX
    for y=1:1:YY
        if M(x,y)==n
            Corx=Corx+x; 
            Cory=Cory+y;
            count=count+1;
        end
    end
    end
    
    Result(n,3)=Corx/count; % x coordinate for the domain
    Result(n,4)=Cory/count; % y coordinate for the domain
end

% save the center of mass data
filename2= 'DomainCenterofMassCoordinates.xlsx';
xlswrite (filename2, Result);
