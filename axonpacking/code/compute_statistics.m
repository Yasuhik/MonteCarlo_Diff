function [FVF, FR, MVF, AVF] = compute_statistics( d, gap, pts, side, g_ratio)
% % author : Tom Mingasson
% stat evaluate from the results of the simulation :
%       - FVF : the fiber volume fraction e.g. the axon (=disk) density
%       - MVF : the myelin volume fraction
%       - AVF : the axon volume fraction
%       - FR  : the fraction of restricted water

disp(' ')
disp('Metrics are being computed...')
disp(' ')

% size of the mask from which metrics are computed. Its area is Atot = Ls * Ls
Ls = sqrt(sum(pi*(d+gap/2).^2))*(4/5);

% resolution in the mask
resolution = 0.2; % um
disp(['pixel size : ', num2str(resolution), ' micro meters'])

% resulting mask size
masksize = ceil(side/resolution);

% FVF mask
FVF_mask = false(masksize);
t = 0:.1:2*pi+0.1;
for id=1:length(d)
    Xfibers = d(id)*cos(t) + pts(1,id);
    Yfibers = d(id)*sin(t) + pts(2,id);
    FVF_mask = FVF_mask | poly2mask(Xfibers/side*masksize, Yfibers/side*masksize, masksize, masksize);
end

% AVF mask
AVF_mask = false(masksize);
for id=1:length(d)
    Xaxons = g_ratio(id)*d(id)*cos(t) + pts(1,id);
    Yaxons = g_ratio(id)*d(id)*sin(t) + pts(2,id);
    AVF_mask = AVF_mask | poly2mask(Xaxons/side*masksize, Yaxons/side*masksize, masksize, masksize);
end

% size of the mask from which metrics are computed. Its area is Atot = Ls * Ls
Ls = sqrt(sum(pi*(d+gap/2).^2))*(4/5)/side*masksize;
Atot = Ls^2;

% center of the mask
Xmin = round(mean(pts(1,:))/side*masksize - Ls/2);
Xmax = round(mean(pts(1,:))/side*masksize + Ls/2);
Ymin = round(mean(pts(2,:))/side*masksize - Ls/2);
Ymax = round(mean(pts(2,:))/side*masksize + Ls/2);

% labeled masks
FVF_mask_trunc = FVF_mask(Xmin:Xmax,Ymin:Ymax);
AVF_mask_trunc = AVF_mask(Xmin:Xmax,Ymin:Ymax);

figure(1000); clf; 
imagesc(0.6.*~FVF_mask_trunc + ~(FVF_mask_trunc - AVF_mask_trunc))
colormap(gray);
axis off;

% labeled areas 
area_fiber = sum(FVF_mask_trunc(:));
area_extra_axonal = Ls*Ls - area_fiber;
area_intra_axonal = sum(AVF_mask_trunc(:));

% compute metrics
FVF = area_fiber / Atot;
FR = area_intra_axonal / (area_intra_axonal + area_extra_axonal); 
AVF = area_intra_axonal / Atot; 
MVF = (area_fiber - area_intra_axonal) / Atot; 

end


