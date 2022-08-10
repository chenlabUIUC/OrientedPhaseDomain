% NAME:     MapStrainGradientVector
% PURPOSE:  This code maps the strain gradient vectors in the cathode nanoparticle
% INPUT:
%           Strain matrix: 'strain_(g002).xlsx'
%           Mask of the cathode particle: 'mask.tif'
% Note:     For generating the mask of the cathode particle, please refer to 
%           the particle morphology in the 4D-STEM data. One can manually delineate the
%           particle shape and make a mask using, for example, imageJ. An example of 
%           mask "mask.tif" is provided in the folder.
% OUTPUT:
%           Strain gradient vector map
% HISTORY:  written by Lehan Yao and Wenxiang Chen, 2022

% Load strain matrix data and strain mask
filename = 'strain_(g002).xlsx';
B = xlsread(filename);
D = imread('mask.tif');
DD=D;
[XX,YY] = size(B);

% Set the background outside the nanoparticle as -5
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B(x,y)=-5;   %background byeond nanoparticle is set as -20
        end
    end
end

% For a better view the strain gradient vector, we choose to bin the strain
% matrix using the imresize function. Note the original strain matrix size
% is [75, 70].
B=imresize(B,[38,35]);
D=logical(imresize(D,[38,35]));
[px,py] = gradient(B);

% For the background outside the nanoparticle, set the strain gradient
% magnitude as -5 and gradient vectors as 0, respectively.
[XX,YY] = size(B);
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
        if x==1 | y==1 | x==XX | y==YY
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% To avoid the sharp gradients at the boundary of the particle, set the
% strain gradient in the pixels close to the particle boundary (within 1 pixel) as 0.
XX=XX-1;
YY=YY-1;
for x=2:1:XX
    for y=2:1:YY
        if D(x-1,y)==0 | D(x+1,y)==0 | D(x,y-1)==0 | D(x,y+1)==0
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% To avoid the sharp gradients at the boundary of the particle, set the
% strain gradient in the pixels close to the particle boundary (within 2 pixels) as 0.
XX=XX-1;
YY=YY-1;
for x=3:1:XX
    for y=3:1:YY
        if D(x-2,y)==0 | D(x+2,y)==0 | D(x,y-2)==0 | D(x,y+2)==0
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% Plot the strain gradient vector in the cathode nanoparticle.
% We first find a max value in strain gradient magnitude, and define a max gradient
% magnitude value as "max_p".
% We then define a scaleing factor of p_s, which is for normalizing the length of
% the gradient arrows. 
% Color scale is defined based on the magnitude of the gradient vectors.
% Gradient vectors are plotted one by one using the quiver function. The
% length of the vectors can be tuned by changing the coefficient of "2.5".
p = sqrt((px).^2+(py).^2);
max_p = 13;
p_s = p ./ ((p/max_p).^.8);
p = 1+floor((p/max_p)*999);
colors = hot(1000);
figure(12);clf;hold on
for xx=1:size(B,1)
    for yy=1:size(B,2)
        quiver(yy,-xx,-py(xx,yy)./p_s(xx,yy).*2.5,px(xx,yy)./p_s(xx,yy).*2.5,'AutoScale','off','color',colors(p(xx,yy),:),'MaxHeadSize',50000)
        %quiver(xx,yy,1,1,'AutoScaleFactor',.9)
        %drawnow
    end
end
axis equal

