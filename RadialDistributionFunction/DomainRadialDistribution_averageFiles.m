% NAME:  DomainRadialDistribution_averageFiles.m 
% PURPOSE:  This script is designed to calculate thge averaged radial distribution of
%                      domains over different MnO2 nanoparticle.
% INPUT:
%           Position matrix: 'Domain center of mass coordinates.xlsx'
%           L: double. Maximum distance.
%           interval: double. Displacement interval.
%           PBC: Switch for periodic boundary condition: "1" for on, "0" for off
% OUTPUT:
%           Plot of pair distribution function for same type of domains and different types of domains
% REFERENCE: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% HISTORY:  written by Chang Qian, 2021/12/29
%

function main()
% Read center of mass of domains.
% Calculate domain radial distribution
% Average over different particles

% Load data
File_list = dir('./');
xls_list = {};
for i = 1:length(File_list)
    if strfind(File_list(i).name,'.xlsx') & strfind(File_list(i).name,'Domain')
        xls_list{end+1} = File_list(i).name;
    end
end
num_files = length(xls_list);

%%%%%%%%%%%% Adjusting parameters %%%%%%%%%%%%%
L = 20; % Maximum displacement, pixel
interval = 1; % Displacement interval, pixel
PBC = 1; % Use periodic boundary condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = struct;
params.L = L; params.interval = interval;
params.PBC = PBC;

RDF_self_final = [];
RDF_diff_final = [];
count_self_total = 0;
count_diff_total = 0;
for i = 1:num_files
    
    pos = xlsread(xls_list{i});
    pos1 = pos(:,1:2); pos1 = pos1(~isnan(pos1(:,1)),:);
    pos2 = pos(:,4:5); pos2 = pos2(~isnan(pos2(:,1)),:);
    [RDF_self, count_self, RDF_diff, count_diff,x]...
        = calcRDF(pos1,pos2,params);
    if isempty(RDF_self_final)
        RDF_self_final = RDF_self;
        RDF_diff_final = RDF_diff;
    else
        RDF_self_final = RDF_self_final + RDF_self;
        RDF_diff_final = RDF_diff_final + RDF_diff;
    end
    count_self_total = count_self_total + count_self;
    count_diff_total = count_diff_total + count_diff;
end

figure(1); clf; hold on;
plot(x,RDF_self_final/count_self_total,'o-','DisplayName','Same type')
plot(x,RDF_diff_final/count_diff_total,'o-','DisplayName','Different type')
xlabel('Displacement (px)')
ylabel('g(r)')
refline(0,1)
legend({'Same type','Different type'})
Outputdata=[];
Outputdata(:,1)=x(:,1);
Outputdata(:,2)=RDF_self_final/count_self_total;
Outputdata(:,3)=RDF_diff_final/count_diff_total;
xlswrite ('pairdistributuionfunction(pixelsize2nm).xlsx', Outputdata);

end


function [RDF_self, count_self, RDF_diff, count_diff,x]...
    = calcRDF(pos1,pos2,params)
    
    L = params.L; interval = params.interval; PBC = params.PBC;

    % Area of particles, for density computation
    Lx1 = max([max(pos1(:,1)),max(pos2(:,1))]);
    Lx2 = min([min(pos1(:,1)),min(pos2(:,1))]);
    Lx = Lx1-Lx2+5;
    Ly1 = max([max(pos1(:,2)),max(pos2(:,2))]);
    Ly2 = min([min(pos1(:,2)),min(pos2(:,2))]);
    Ly = Ly1-Ly2+5;
    % Lx = max([max(pos1(:,1)),max(pos2(:,1))]);
    % Ly = max([max(pos1(:,2)),max(pos2(:,2))]);
    
    rdf1 = calcRadialDistPBC2DSelect(pos1, pos1, Lx, Ly, L, interval, PBC);
    rdf2 = calcRadialDistPBC2DSelect(pos1, pos2, Lx, Ly, L, interval, PBC);
    rdf3 = calcRadialDistPBC2DSelect(pos2, pos1, Lx, Ly, L, interval, PBC);
    rdf4 = calcRadialDistPBC2DSelect(pos2, pos2, Lx, Ly, L, interval, PBC);
    
    x = rdf1(:,1);
    RDF_self = rdf1(:,2) + rdf4(:,2);
    RDF_diff = rdf2(:,2) + rdf3(:,2);
    count_self = size(pos1,1)+size(pos2,1);
    count_diff = size(pos1,1)+size(pos2,1);
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