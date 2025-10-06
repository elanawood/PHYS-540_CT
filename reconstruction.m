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

% 2. consider 3 views, 0 degrees (top-bottom), 45 degrees (diagonal), and 90
% degrees (left-right)
%   - calculate the projection for each view 
%   - plot each projection (this is a sinogram)

%define empty arrays for projection values
view0 = zeros(1,21);
view90 =zeros(1,21);
n=size(medium,1);
view45=zeros(1,2*n-1);

%define max pixel lengths xi for each view
x0=1;
x90=1;
x45=sqrt(2);

% calculate projections for 0 and 90 degrees
for col=1:21
    for row=1:21
        view0(col)= view0(col)+medium(row,col)*x0;
        view90(col)=view90(col)+medium(col,row)*x90;
    end
end

% calculate projections for 45 degrees
% flip the matrix vertically so that the antidiagonals of 'medium' become
% the main diagonals of 'B'

B=flipud(medium);

for k=-(n-1):(n-1)
    idx=k+n;
    view45(idx)=sum(diag(B,k))*x45;
end




