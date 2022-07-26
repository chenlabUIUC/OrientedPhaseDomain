% NAME:     CalculateDomainStatistics
% PURPOSE:  This code calculates the statistics of oriented phase domains in the particle 
% INPUT:
%           Domain map: 'DomainMap.xlsx'
%           Mask: 'mask.tif'
% OUTPUT:
%           An excel file containing statistics of the oriented phase domains:
%           'DomainStatistics.xlsx' (total number of domains, particle area, total area of domains)
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019

% Load data the domain distribution map
filename = 'DomainMap.xlsx';
B = xlsread(filename);
C=B;
D=B;
[XX,YY] = size(B);

% Construct matrixes for the [100]t and [111]t domains, respectively
for x=1:1:XX
    for y=1:1:YY
        if B(x,y)~=2
            C(x,y)=0; % select oriented phase domain of [100]t orientation
        end
        if B(x,y)~=1
            D(x,y)=0; % select oriented phase domain of [111]t orientation
        end
    end
end

% Label connected components in 2-D binary image
L = bwlabel(C,4);
M = bwlabel(D,4);

% Optional: show the labelled oriented phase domains
%{
figure
P=imagesc(L);
hold on
figure
P=imagesc(M);
hold off
%}

% Apply frequency count function to calculate the area of each oriented phase domain
LL=L(:);
MM=M(:);
LLmax=max(LL);
MMmax=max(MM);
[x1,y1]=hist(LL,2*LLmax);  % x1 returns the list of area (unit: pixel) of each [110]t domain; It contains the background and zeros
[x2,y2]=hist(MM,2*MMmax);  % x2 returns the list of area (unit: pixel) of each [110]t domain; It contains the background and zeros

% calculate the area of [100]t oriented phase domain
DomainAreaA = x1*4; % x1*4 convert the unit from pixel^2 to nm^2 (1 pixel^2=4 nm^2)
CountA = numel(DomainAreaA);
DomainAreaA2 = DomainAreaA(2:CountA); % remove the background, which is usually the first element in the array
DomainAreaA2 = nonzeros(DomainAreaA2); % remove the zeros generated in the frequency count
% calculate the area of [111]t oriented phase domain
DomainAreaB = x2*4; % x1*4 convert the unit from pixel^2 to nm^2 (1 pixel^2=4 nm^2)
CountB = numel(DomainAreaB);
DomainAreaB2 = DomainAreaB(2:CountB); % remove the background, which is usually the first element in the array
DomainAreaB2 = nonzeros(DomainAreaB2); % remove the zeros generated in the frequency count

summary=[DomainAreaA2;DomainAreaB2];

NumberDomain = numel(summary);      % Total domain number in the particle
TotalDomainArea = sum(summary);  % Total domain area (nm^2) in the particle

mask = imread('mask.tif');
ParticleArea = 4*sum(mask(:) == 255); % Particle area (nm^2)

% Output the Number of tetragonal domains, the cathode particle area, and
% the total area of the tetragonal domains in the cathode particle
result1 = {'NumberOfDomains', 'ParticleArea(nm2)', 'TotalDomainArea(nm2)'};
result2 = [NumberDomain, ParticleArea, TotalDomainArea;] ;
result = [result1;num2cell(result2)];
xlswrite('DomainStatistics.xlsx',result);
