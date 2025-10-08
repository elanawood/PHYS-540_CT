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
% the main diagonals of 'B', working bottom left -> top right in final
% result

B=flipud(medium);

for k=-(n-1):(n-1) % map k-> 1..(2n-1), sum along main diagonal offset k
    idx=k+n;
    view45(idx)=sum(diag(B,k))*x45;
end

% plot the projections - stack them into a small image
% 0 degree projection
figure;
imagesc(view0);axis image; colormap("gray");colorbar;
title("Projection at 0 degrees");
xlabel("Detector value");
set(gca,'Ytick',[]);

% 90 degree projection
figure;
imagesc(view90(:));axis image; colormap("gray");colorbar;
title("Projection at 90 degrees");
ylabel("Detector value");
set(gca,'Xtick',[]);


%45 degree projection
figure;
imagesc(view45);axis image; colormap("gray");colorbar;
title("Projection at 45 degrees, full length");
xlabel("Detector value");
set(gca,'Ytick',[]);

% pad view45 to length 21 to plot sinogram with no backprojection
view45_padded=zeros(1,21);
middle=n;
s_idx=middle-floor(n/2);
view45_padded(:)=view45(s_idx:s_idx+n-1);

sinogram_padded=[view0;view45_padded;view90];

figure;
imagesc(sinogram_padded);
colormap("gray");colorbar;
yticks(1:3);yticklabels({'0 degrees','45 degrees (centered to 21 bins)','90 degrees'});
xlabel('Detector Value'); ylabel('Angle');
title("Basic Sinogram (3 angles), No Backprojection");


% backproject 0 and 90 degree view vectors across rows and columns
backprojection0=repmat(view0,n,1);
backprojection90=repmat(view90(:),1,n);

% backproject 45 degrees, using the padded 45 degree view with 21 elements
% from bottom left -> top right

backprojection45 = zeros(size(medium));
start_index=(n+1)/2;

for k=1:n
    t=start_index+ (k-1);

    % anti diagonal pixels
    i_min=max(1,t+1-n);
    i_max=min(n,t);

    val=view45_padded(k); % the same values appear across the diagonal

    % iterate bottom left to top right

    for i=i_max:-1:i_min
        j=t+1-i;
        backprojection45(i,j)=backprojection45(i,j)+val;
    end
end

%sum all backprojections into a single image
backprojection= backprojection0+backprojection45+backprojection90;
% Display the final backprojection result
figure;
imagesc(backprojection);
colormap("gray");
colorbar;
title("Simple Backprojection");
xlabel("X-axis (pixels)");
ylabel("Y-axis (pixels)");
axis image;


% filtered backprojection
% define filter:
F=[-0.1074,0.1368,-0.3398,0.6000,-0.3398,0.1368,-0.1074];
% Convolve each projection vector with the filter F
filteredView0 = conv(view0, F, 'same');
filteredView90 = conv(view90, F, 'same');
filteredView45 = conv(view45_padded, F, 'same');


% backproject the filtered projection vectors
filteredbackprojection0=repmat(filteredView0,n,1);
filteredbackprojection90=repmat(filteredView90(:),1,n);

filteredbackprojection45 = zeros(size(medium));
start_index=(n+1)/2;

for k=1:n
    t=start_index+ (k-1);

    % anti diagonal pixels
    i_min=max(1,t+1-n);
    i_max=min(n,t);

    val=filteredView45(k); % the same values appear across the diagonal

    % iterate bottom left to top right

    for i=i_max:-1:i_min
        j=t+1-i;
        filteredbackprojection45(i,j)=filteredbackprojection45(i,j)+val;
    end
end

% construct the filtered backprojection image
filteredbackprojection = filteredbackprojection0 + filteredbackprojection45 + filteredbackprojection90;

% Display the filtered backprojection result
figure;
imagesc(filteredbackprojection);
colormap("gray");
colorbar;
title("Filtered Backprojection");
xlabel("X-axis (pixels)");
ylabel("Y-axis (pixels)");
axis image;

% calculating a 1/r blurring for each image
% center index (11,11)
cx=11;cy=11;

% radius map in pixels
[xg,yg]= meshgrid(1:n,1:n);
R=sqrt((xg-cx).^2+(yg-cy).^2);

% bin by integer radius, rings are 1 pixel wide
rbin=floor(R);
rmax=max(rbin(:));

% helper function to calculate the ring-mean
ring_mean=@(img) arrayfun(@(k) mean(img(rbin==k)),0:rmax);

%calculate the radial mean for each image
rMean_bpg=ring_mean(backprojection);
rMean_fbpg=ring_mean(filteredbackprojection);



% display results graphically
figure;
plot(r, rMean_bpg,  'o-', 'LineWidth', 1.5, 'DisplayName', 'Simple Backprojection'); hold on;
plot(r, rMean_fbpg, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Filtered Backprojection');

grid on;
xlabel('Ring radius (pixels)');
ylabel('Average pixel value in ring');
title('Radial Blur Comparison (Unnormalized)');
legend('Location','northeast');

% display results as 2D grayscale images
ringImage_simple  = zeros(n,n);
ringImage_filtered = zeros(n,n);
ringImage_reference = zeros(n,n);

for k = 0:rmax
    ringImage_simple(rbin == k)   = rMean_bpg(k+1);  
    ringImage_filtered(rbin == k) = rMean_fbpg(k+1);
    ringImage_reference(rbin == k)= one_over_r(k+1);
end

figure;
subplot(1,3,1); imagesc(ringImage_simple); colormap(gray); axis image off;
title('Simple BP — Radial Blur Map');

subplot(1,3,2); imagesc(ringImage_filtered); colormap(gray); axis image off;
title('Filtered BP — Radial Blur Map');

subplot(1,3,3); imagesc(ringImage_reference); colormap(gray); axis image off;
title('1/r Reference Map');