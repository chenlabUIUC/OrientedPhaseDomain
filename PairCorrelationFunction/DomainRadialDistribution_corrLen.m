% NAME:  DomainRadialDistribution_corrLen.m
% PURPOSE:  This script is designed to calculate radial distribution of
%                      domains in MnO2 cathode material.
% INPUT:
%           Position matrix: 'Domain center of mass coordinates.xlsx'
%           L: double. Maximum distance.
%           interval: double. Displacement interval.
%           PBC: Switch for periodic boundary condition: "1" for on, "0" for off
% OUTPUT:
%           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% REFERENCE: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% HISTORY:  written by Chang Qian, 2021/12/29
%

function main()
% Read center of mass of domains.
% Calculate domain radial distribution

% Load data
pos = xlsread('Domain center of mass coordinates.xlsx');
pos1 = pos(:,1:2); pos1 = pos1(~isnan(pos1(:,1)),:);
pos2 = pos(:,4:5); pos2 = pos2(~isnan(pos2(:,1)),:);
pos_mat = [pos1,ones(size(pos1,1),1);pos2,ones(size(pos2,1),1)*2]

%%%%%%%%%%%% Random data%%%%%%%%%%%%%%%
% n_type1 = 80; % number of type1 data points
% n_type2 = 40; % number of type 2 data points
% area_max = 100; % data expand in this area.
% pos1 = rand(n_type1,2)*area_max;
% pos2 = rand(n_type2,2)*area_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Adjusting parameters %%%%%%%%%%
L = 20; % Maximum displacement, pixel
interval = 1; % Displacement interval, pixel
PBC = 1;  %periodic boundary condition: "1" means on, "0" means off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Area of particles, for density computation
Lx1 = max(pos1(:,1))-min(pos1(:,1)); Lx2 = max(pos2(:,1))-min(pos2(:,1)); Lx = max(Lx1,Lx2) + 10; 
Ly1 = max(pos1(:,2))-min(pos1(:,2)); Ly2 = max(pos2(:,2))-min(pos2(:,2)); Ly = max(Ly1,Ly2) + 10;  %% Might need to adjust the value 10 %%
% Lx = area_max; Ly = area_max;

% Using periodic boundary condition?
% PBC = 1; 
% Lx1 = max(pos1(:,1)); Lx2 = max(pos2(:,1)); Lx = max(Lx1,Lx2); 
% Ly1 = max(pos1(:,2)); Ly2 = max(pos2(:,2)); Ly = max(Ly1,Ly2);


figure(1); clf; hold on; box on; axis equal;
set(gca,'Ydir','reverse');
scatter(pos1(:,2),pos1(:,1),'ro','filled')
scatter(pos2(:,2),pos2(:,1),'bo','filled')

figure(2); clf; hold on;
% rdf = calcRadialDistPBC2DSelect(pos1, pos1, Lx, Ly, L, interval, PBC);
out_mat = calcCorrLenPBC2DSelect(pos_mat, Lx, Ly, L, interval, PBC);
CorrLen = out_mat.mat;
CorrLen(:,4) = CorrLen(:,2)./CorrLen(:,3);
x = CorrLen(:,1);
y = CorrLen(:,4);
figure(2); clf; hold on;
plot(x,y,'o-','LineWidth',2)
refline(0,0)


%%%%%%%%%%%% Direct plotting %%%%%%%%%%
% plot(rdf(:,1),rdf(:,2)/size(pos1,1),'o-','DisplayName','type1-type1')
% data=[];
% data(:,1)=rdf(:,1);
% DomainNumber1=size(pos1,1);
% data(:,2)=rdf(:,2)/size(pos1,1);
% rdf = calcRadialDistPBC2DSelect(pos1, pos2, Lx, Ly, L, interval, PBC);
% plot(rdf(:,1),rdf(:,2)/size(pos1,1),'o-','DisplayName','type1-type2')
% data(:,3)=rdf(:,2)/size(pos1,1);
% xlabel('Displacement (px)')
% ylabel('g(r)')
% legend
% 
% figure(3); clf; hold on;
% rdf = calcRadialDistPBC2DSelect(pos2, pos1, Lx, Ly, L, interval, PBC);
% plot(rdf(:,1),rdf(:,2)/size(pos2,1),'o-','DisplayName','type2-type1')
% data(:,4)=rdf(:,2)/size(pos2,1);
% DomainNumber2=size(pos2,1);
% rdf = calcRadialDistPBC2DSelect(pos2, pos2, Lx, Ly, L, interval, PBC);
% plot(rdf(:,1),rdf(:,2)/size(pos2,1),'o-','DisplayName','type2-type2')
% data(:,5)=rdf(:,2)/size(pos2,1);
% xlabel('Displacement (px)')
% ylabel('g(r)')
% legend
% 
% data(:,6)=(data(:,2)*DomainNumber1+data(:,5)*DomainNumber2)/(DomainNumber1+DomainNumber2);
% data(:,7)=(data(:,3)*DomainNumber1+data(:,4)*DomainNumber2)/(DomainNumber1+DomainNumber2);
% %xlswrite ('data.xlsx', data);
% 
% figure(4); clf; hold on;
% plot((data(:,1)+0.5)*2,data(:,6),'o-','DisplayName','Same domain')
% plot((data(:,1)+0.5)*2,data(:,7),'o-','DisplayName','Diff domain')
% xlabel('Displacement (nm)')
% ylabel('g(r)')
% legend

%%%%%%%%%%%% Data output %%%%%%%%%%
% 1st column as distance (pixel), 
% 2nd column as sum of +1 same orientation or -1 different orientation, 
% 3rd column as count of domains at the distance
% 4th column as the ratio 2nd/3rd column.
xlswrite ('CorrLen', CorrLen);   
end


function rdf_final = calcRadialDistPBC2DSelect( pos_ref, pos, Lx, Ly, L, interval, PBC)
% NAME:     calcRadialDistPBC2D.m
% PURPOSE:  calculate radial distribution function in two dimension system.
%           Use periodic boundary conditions to fix the finite size effect
%           in the system. 
%           The system should be in a rectangular shape for it uses
%           periodic boundary condition to calculate the distance between
%           two particles.
%           The particles should be able to move freely in the whole frame
%           for the ideal gas density is defined as the number of particles
%           over the system size.
%           The range of radius calculated is [0, 0.3*(Lx^2 + Ly^2)^(1/2)].
%           This is kept the same with the referred code.
% http://www.mathworks.com/matlabcentral/fileexchange/46576-radialdistribution2d-m
% http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
% INCLUDE:  displacementPBC = distPBC2D(displacement2D, Lx, Ly)
%               this function is only used in this calculation so it is
%               attached at the last. This will return the distance between 
%               two particles under periodic boundary condition.
%           h = addHist(h,data)
%               this function is only used in this calculaiton so it is
%               attached at the last. This function will add one element
%               towards the existing histogram.
% INPUT:    pos: particle position list. N * 2 matrix, N : total number of
%               particles. The row number is the particle index.
%               pos(:,1): column array of x position.
%               pos(:,2): column array of y position.
%           Lx, Ly: system length in x and y direction.
%           n_interval: the number of intervals to divide in radial
%               distribution calculation.
%           list
% OUTPUT:   rdf: radial distribution function. n_interval * 2 matrix.
%               rdf(:,1): radius in pixel.
%               rdf(:,2): g(r) values corresponding to specific radius.
% HISTORY:  written by zihao, 20151105
%           modified by zihao, 20170614, change the calculation to allow
%               for single particle radial distribution function.
%           modified by Chris, 20211229, prepare data matrix for correlation length
%           analysis for Wenxiang
%%
% initialize rdf
rdf = struct;
rdf.count = 0;
rdf.range = [0 L];
rdf.increment = interval;
% calculate rdf based on PBC assumption
% Loop over pairs and determine the distribution of distances
[nPart_ref,~] = size(pos_ref);
[nPart,~] = size(pos);
for partA = 1:nPart_ref
    for partB = 1:nPart
        if isequal(pos_ref(partA,:),pos(partB,:))
            continue;
        end
        % Calculate particle-particle distance                             
        displacement = pos_ref(partA,:) - pos(partB,:); 
        % Account for PBC (assuming 2D)  % Periodic Boundary Correction
        if PBC
            displacement = distPBC2D(displacement,Lx,Ly);
        end
        dx = displacement(1);   dy = displacement(2);
        r = sqrt(dx * dx + dy * dy);
        % Add to g(r) if r is in the right range [0 0.3*L]
        if (r < L)
            rdf = addHist(rdf,r);
%             if r < 10
%                 disp(partA); disp(partB); pause;
%             end
        end
    end
end
% each bin should be normalized according to its volume
nBins = size(rdf.values,2);%second dimension, rdf is a row vector
rho = nPart/(Lx*Ly);
for bin = 1:nBins
%   the area of the ring is from (bin-1) * increment to bin * increment,
%   not r to r + increment. r is the average radius of the ring.
    rVal = rdf.increment * (bin - 1);
    next_rVal = rVal + rdf.increment;
    areaBin = pi*next_rVal^2 - pi*rVal^2; 
    % Calculate the number of particles expected in this bin in
    % the ideal case
    nIdeal = areaBin*rho;
    % Normalize the bin
    rdf.histo(bin) = rdf.histo(bin) / nIdeal;
end
% The radial distribution function should be normalized.
rdf.histo = rdf.histo;
% out put final rdf result.
rdf_final = [rdf.values',rdf.histo']; 
end

function out_mat = calcCorrLenPBC2DSelect( pos_mat, Lx, Ly, L, interval, PBC)
% NAME:     calcCorrLenPBC2D.m
% 
% HISTORY:  written by zihao, 20151105
%           modified by zihao, 20170614, change the calculation to allow
%               for single particle radial distribution function.
%           modified by Chris, 20211229, prepare data matrix for correlation length
%           analysis for Wenxiang
%%
% initialize rdf
out_mat = struct;
out_mat.count = 0;
out_mat.range = [0 L];
out_mat.increment = interval;
out_mat.mat = [interval:interval:L]';

% calculate rdf based on PBC assumption
% Loop over pairs and determine the distribution of distances
[nPart,~] = size(pos_mat);
out_mat.mat = [out_mat.mat,zeros(size(out_mat.mat,1),2)];
for partA = 1:nPart
    for partB = 1:nPart
        if isequal(pos_mat(partA,:),pos_mat(partB,:))
            continue;
        end
        % Calculate particle-particle distance                             
        displacement = pos_mat(partA,:) - pos_mat(partB,:); 
        % Account for PBC (assuming 2D)  % Periodic Boundary Correction
        if PBC
            displacement = distPBC2D(displacement,Lx,Ly);
        end
        dx = displacement(1);   dy = displacement(2);
        r = sqrt(dx * dx + dy * dy);
        % Add to g(r) if r is in the right range [0 0.3*L]
        if (r < L)
%             rdf = addHist(rdf,r);
%             if r < 10
%                 disp(partA); disp(partB); pause;
%             end
            ind = ceil(r/interval);
            out_mat.mat(ind,3) = out_mat.mat(ind,3)+1;
            if pos_mat(partA,3)==pos_mat(partB,3)
                out_mat.mat(ind,2) = out_mat.mat(ind,2)+1;
            else
                out_mat.mat(ind,2) = out_mat.mat(ind,2)-1;
            end
        end
    end
end
end

%%
function d = distPBC2D(d,Lx,Ly)
% NAME: distPBC2D.m
% PURPOSE:  calculate the distance between two particles in periodic
%           boundary condition. The nearest distance between two particles
%           including their images.
% INCLUDE:  none
% INPUT:    d: displacement. (x, y)
%           Lx: system size in x direction.
%           Ly: system size in y direction.
% OUTPUT:   d: displacement in PBC condition. (x, y)
% http://www.mathworks.com/matlabcentral/fileexchange/46575-distpbc2d-m
% HISTORY:  written by zihao, 20151105

%   Calculate the half box size in each direction
    hLx = Lx/2.0;    hLy = Ly/2.0;
%   Calculate the distance of the nearest particle distance   
    if d(1) > hLx
        d(1) = d(1) - Lx;
    elseif d(1) < -hLx
        d(1) = d(1) + Lx;
    end
    if d(2) > hLy
        d(2) = d(2) - Ly;
    elseif d(2) < -hLy
        d(2) = d(2) + Ly;
    end
end
%%
function h = addHist(h,data)
% NAME:     addHist.m
% PURPOSE:  this program will add an element to an existing histogram. If
%           the element is the first elecment adding towards the histogram,
%           it will initialize the histomgram.
% http://www.mathworks.com/matlabcentral/fileexchange/46576-radialdistribution2d-m
% INCLUDE:  none.
% INPUT:    h:  a struct used to store information for the histogram.
%           data:   the element adding towards the histogram.
% OUTPUT:   h:  the struct of histogram with the data adding into it.
%               h.values: the center value of each column.
%               h.histo: the number of data in each column.
% HISTORY:  written by zihao, 20151106
    if (h.count == 0)
        % Determine the number of bins by evaluating the histogram''s range and increment size  
        nBins = ceil((h.range(2)-h.range(1))/h.increment);
        % Adjust the histogram''s range
        % Useful if the total range is not an exact multiple of the increment size
        h.range(2) = h.range(1) + nBins * h.increment;
        % Set all bins to zero
        h.histo = zeros(1,nBins);
        % Set the values vector to be in the center of each bin
        h.values = 1:nBins;
        h.values = h.range(1) + h.increment*(h.values-0.5);
    end
%   Check if the data is in the range.
    if (data > h.range(1) && data <= h.range(2)) % Make sure the data fits the range
        % Find the right bin position
        binIndex = ceil((data-h.range(1))/h.increment);
        h.histo(binIndex) = h.histo(binIndex)+1;
        h.count = h.count+1;%total count
    else
        disp(['histogram- Value out of range:',num2str(data)]);
        return
    end
end

