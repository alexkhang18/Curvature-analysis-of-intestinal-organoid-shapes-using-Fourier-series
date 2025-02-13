clc
clear all
close all

%% Fourier series settings

total_number_of_series = 10;

%% image resolution

x_res = 369.28/800; % microns/pixel
y_res = x_res;

%% Load in files of interest

% .tif orgnanoid body image
I = imread('organoid_brightfield.tif');
I = imadjust(I);
I = cat(3,I,I,I);

% read in organoid body segmentation 
organoid = imread('organoid_segmentation.png');

% read in lysozyme segmentation 
lysozyme = imread('lysozyme_segmentation.png');

%% Image registration, Fourier fit, and curvature computation 
    
% extracts boundary point of current shape
BWoutline = organoid;
BWoutline = imclearborder(BWoutline);
BWoutline = BWoutline>0;
BWfill = imfill(BWoutline,'holes');
positive = find(BWfill>0);
SA = length(positive) * x_res * y_res;
BWoutline = bwperim(BWoutline);
[x,y] = find(BWoutline>0);
points = [x,y];

% algorithm to sort points
points2 = points; % copy of original points
points_sort = [x(1),y(1)]; % starts points_sort at x(1),y(1)
points2(1,:) = []; % clears points as they are added to points_sort
for i = 1:length(x)-1 % iterates through remaining points
    distance = vecnorm(points2-points_sort(end,1:2),2,2); % computes distance between all remaining points
    distance(distance==0) = Inf; % prevents current point from being counted again
    indx = find(distance == min(distance)); % finds nearest point
    points_sort = [points_sort;points2(indx(1),:)]; % appends nearest point to points_sort
    points2(indx(1),:) = []; % clears points as they are added to points_sort
end

% converts boundary from pixels to um and from image indexing to cartesian
% coordinates 
image_to_points = [points_sort(:,2)*y_res,(size(BWoutline,1)-points_sort(:,1))*x_res];

% adds the first point to the end of array to ensure a closed loop
points_um = [image_to_points;image_to_points(1,:)];

% finds arc-length of distance around organoid shape
dt = diff(points_um);
dis = zeros(length(dt),1);
for i = 1:length(dt)
    if i == 1
        dis(i) = norm(dt(i,:));
    else
        dis(i) = dis(i-1) + norm(dt(i,:));
    end
end

% prepares arc-length data for Fourier series fit
dis = [0;dis];
dis = normalize(dis,'range')*2*pi;

% concatenates x and y points separately 
points_um_x = points_um(:,1);
points_um_y = points_um(:,2);

% computes Fourier series coefficients
[a] = fourier_series_fit(dis,points_um_x,total_number_of_series);
[b] = fourier_series_fit(dis,points_um_y,total_number_of_series);

% Fourier series evaluated at these theta values
thetas = linspace(0,2*pi,10000);

% Evaluation of Fourier series fit
[x_predicted,fx,fxx] = fourier_series(a,thetas);
[y_predicted,fy,fyy] = fourier_series(b,thetas);

% Evaluation of curvature kappa
kappa = (fx.*fyy - fy.*fxx)./((fx.^2+fy.^2).^(3/2));

% ensures color scheme matches curvature kappa 
col = kappa;

% variable used for plotting later
z_predicted = zeros(size(x_predicted));    

%% arc-length vs kappa 

% finds arc-length of distance around Fourier series fit of organoid shape
points_predicted = [x_predicted,y_predicted];
dt = diff(points_predicted);
dis = zeros(length(dt),1);
for i = 1:length(dt)
    if i == 1
        dis(i) = norm(dt(i,:));
    else
        dis(i) = dis(i-1) + norm(dt(i,:));
    end
end
dis = [0;dis];

% curvature bins
k_intervals = -0.305:0.01:0.305;

% iterate through all the kappa intervals to find segments of the organoid
% boundary that falls into the bin
for kk = 1:length(k_intervals)-1

    % find kappa values which fall into the current range
    indx_temp = find(kappa >= k_intervals(kk) & kappa <= k_intervals(kk+1));

    % finds the arc-length values for the identified kappa values
    dis_temp = dis(indx_temp);

    % finds the difference in arc-length between identified kappa values 
    diff_dis_temp = diff(dis_temp);

    % removes sections with unreasonably large distances which signify the 
    % distance between two non-neighboring points along the organoid 
    % boundary which fall into the current curvature range 
    diff_dis_temp(diff_dis_temp>max(diff(dis))) = [];

    % removes redundant sections that are likely points that were counted
    % twice (very unlikely this code segment is needed)
    diff_dis_temp(diff_dis_temp<min(diff(dis))) = [];

    %computes middle bin value for plotting
    k_mid(kk) = mean([k_intervals(kk),k_intervals(kk+1)]);

    % sums arc-length along organoid body with curvature that falls into
    % the current bin
    arc_length(kk) = sum(diff_dis_temp);

end
 
% concatenates results for plotting
k_vs_arclength = [k_mid',arc_length'];

% arc length vs curvature
figure
bar(k_mid,arc_length);
xlim([-0.2 0.2])
ylim([0 160])
xlabel('$\kappa$ ($\mu$m)','Interpreter','Latex')
ylabel('arc length ($\mu$m)','Interpreter','Latex')
hold off;
exportgraphics(gcf,'arclength_vs_kmid.png','Resolution',300)

%% distance-based weighted average interpolation of curvature from boundary to inside the organoid body

% converts boundary from cartesian coordinates back to pixels and image 
% indices
x_predicted = x_predicted * (1/x_res);
y_predicted = y_predicted * (1/x_res);
y_predicted = size(BWfill,1) - y_predicted;

% interpolates kappa onto image of organoid boundary 
points = [x_predicted,y_predicted]; % concatenates points
TEST = NaN(size(BWfill)); % pre-allocation of curvature image
x_glob = round(points(:,1)); % x value of boundary pixels
y_glob = round(points(:,2)); % y value of boundary pixels

% iterate through boundary points and fill in curvature image and creates 
% kappa global variable for interpolation scheme
for j = 1:length(x_glob) 
    TEST(y_glob(j),x_glob(j)) = kappa(j);
    kappa_glob(j,:) =  kappa(j);
end

% creates mask of organoid body that has been eroded by two pixels to avoid
% edge effects in calculations 
BW = TEST>-1000; % all curvature values are above -1000 so this ensures the entire boundary is accounted for
BW = imfill(BW,'holes'); % fills holes to get organoid body
BWperim = bwperim(BW); % finds perimeter values
BW(BWperim) = 0; % erodes perimeter by setting values to 0
BWperim = bwperim(BW); % finds perimeter values
BW(BWperim) = 0; % erodes perimeter by setting values to 0
lysozyme(~BW) = 0; % erases all lysozyme signal outside of organoid body as it is likely not true signal

% finds pixel locations of all pixels inside organoid body 
[row,column] = find(BW > 0);

% pixel locations of boundary points 
bound = [y_glob,x_glob];

% pre-allocation of curvature range and color scale
kappa_range = linspace(-0.1,0.1,2^16);
c = jet(2^16);

% pre-allocation of an image based on brightfield image
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);

% iterates through each pixel inside organoid body
for j = 1:length(row)

    % computes distance between current pixel and all boundary pixels
    distance = vecnorm(bound - [row(j),column(j)],2,2);

    % computes weight which is inversely related to distance squared
    weight = (1./(distance.^2));

    % interpolated curvature value for current pixel based on weighted
    % average of curvature values at the boundary pionts
    value = sum((weight.*kappa_glob))/sum(weight);

    % linear interpolation of color scheme to find color corresponding to
    % the kappa value computed for the current pixel
    vq = interp1(kappa_range,c,value);
    
    % fills in TEST variable which is the interpolated curvature image
    TEST(row(j),column(j)) = value;

    % fills in brightfield image with interpolated curvature mask. converts
    % the interpolated color values to 16-bit image using 2^16
    redChannel(row(j),column(j)) = round(2^16*vq(1));
    greenChannel(row(j),column(j)) = round(2^16*vq(2));
    blueChannel(row(j),column(j)) = round(2^16*vq(3));

end

% concatenates final image
rgbImage = cat(3, redChannel, greenChannel, blueChannel);

% curvature interpolation into organoid body
figure 
imshow(rgbImage); hold on;
colormap('jet')
clim([-0.10 0.10]); hold off;
title('interpolated boundary curvature into organoid body')
hold off;
exportgraphics(gcf,'interpolated_curvature_organoid.png','Resolution',300)

%% Applying interpolated curvature values to regions with positive Lysozyme staining

% finds pixel locations with positive lysozyme staining
[row,column] = find(lysozyme>0);

% pre-allocation of curvature range and color scale
kappa_range = linspace(-0.1,0.1,2^16);
c = jet(2^16);

% pre-allocation of an image based on brightfield image
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);

% pre-allocation 
kappa_temp = [];

% loop through every positive lysozyme pixel
for i = 1:length(row)

    % records the curvature values at each positive lysozyme pixel
    kappa_temp(i) = TEST(row(i),column(i));

    % linear interpolation of color scheme to find color corresponding to
    % the kappa value computed for the current pixel
    vq = interp1(kappa_range,c,kappa_temp(i));

    % fills in brightfield image with interpolated lysozyme curvature mask. converts
    % the interpolated color values to 16-bit image using 2^16
    redChannel(row(i),column(i)) = round(2^16*vq(1));
    greenChannel(row(i),column(i)) = round(2^16*vq(2));
    blueChannel(row(i),column(i)) = round(2^16*vq(3));

end

% concatenates final image
rgbImage = cat(3, redChannel, greenChannel, blueChannel);

% interpolated curvature into lysozyme segmentation 
figure
imshow(rgbImage); hold on;
colormap('jet')
clim([-0.10 0.10]);
title('interpolated curvature in positive lysozyme pixels')
hold off;
exportgraphics(gcf,'interpolated_curvature_lysozyme_pixels.png','Resolution',300)

%% plots for visualization

pause(1)

% fourier series fit and lysozyme segmentation overlayed on top of organoid body, colored by
% curvature
figure
imshow(labeloverlay(I,lysozyme,'Transparency',0.5,'Colormap',[1,0,1])); hold on;
surface([x_predicted';x_predicted'],[y_predicted';y_predicted'],[z_predicted';z_predicted'],[col';col'],...
        'facecol','no',...
        'edgecol','flat',...
        'linew',10); 
colormap('jet')
clim([-0.10 0.10]); hold off;
title('boundary curvature and positive lysozyme pixels')
hold off;
exportgraphics(gcf,'boundary_curvature_and_lysozyme_pixels.png','Resolution',300)

% probability density estimate of interpolated curvature values in lysozyme positive
% regions
k_range = linspace(-0.2,0.2,1000);
[f,xi] = ksdensity(kappa_temp,k_range); 
figure
plot(xi,f);
xlim([-0.1 0.1])
xlabel('$\kappa$ ($\mu$m)','Interpreter','Latex')
ylabel('probability density estimate')
title('curvature values associated with lysozyme pixels')
hold off;
exportgraphics(gcf,'PDE_lysozyme_curvatures.png','Resolution',300)

% colorbar
figure        
c = colorbar;
c.Label.String = '\kappa (\mum^-1)';
colormap('jet')
clim([-0.10 0.10]); hold off;
set(gca,'Visible','off')
exportgraphics(gcf,'colorbar.png','Resolution',300)
