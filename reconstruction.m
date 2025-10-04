% Script to run simple and filtered reconstruction on a test object

% Steps used
% 1. Construct test pattern - non attenuating medium with an attenuating
% block in the center. 
%   - Medium is a square matrix, 21 by 21 pixels in size.
%   - Medium has linear attenuation coeff u=0 mm^-1 
%   - Central 3x3 area is attenuating u=5 mm^-1
%   - Each pixel is 1x1 mm

% define medium
medium=zeros(21,21); % Define 21x21 matrix to represent medium
medium(10:12, 10:12) = 5; % Set the central 3x3 area to the attenuating value