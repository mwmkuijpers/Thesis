close all 
clear all
warning off
i = num2str(1);
imgBottom = imgaussfilt(imread(['C:\Users\Marly\Dropbox\ThesisMatlab\BottomCamera\img' i '.jpg']),10);
imgFront = (imread(['C:\Users\Marly\Dropbox\ThesisMatlab\FrontCamera\img' i '.jpg']));

% separate channels
redChannelBottom = double(imgBottom(:,:,1));
greenChannelBottom = double(imgBottom(:,:,2));
blueChannelBottom = double(imgBottom(:,:,3));

redChannelFront = double(imgFront(:,:,1));
greenChannelFront = double(imgFront(:,:,2));
blueChannelFront = double(imgFront(:,:,3));

% mean of each channel
redMeanBottom = mean(redChannelBottom(:));
greenMeanBottom = mean(greenChannelBottom(:));
blueMeanBottom = mean(blueChannelBottom(:));

% Determine Euclidian distance for every pixel
% The Euclidean distance between points s and t is the length of the line segment connecting them:
% d_(st)^2 = (x_s-x_t)*(x_s-x_t)' 
% or: d_st = sqrt((s1-t1)^2+(s2-t2)^2...)
% determine color difference 
redB = redMeanBottom*ones(size(redChannelFront,1),size(redChannelFront,2));
greenB = greenMeanBottom*ones(size(greenChannelFront,1),size(greenChannelFront,2));
blueB = blueMeanBottom*ones(size(blueChannelFront,1),size(blueChannelFront,2));
redDistEucl = sqrt((redChannelFront-redB).^2);
greenDistEucl = sqrt((greenChannelFront-greenB).^2);
blueDistEucl = sqrt((blueChannelFront-blueB).^2); 

% set thresholds

redThreshold = 20;%mean(redDistEucl(:));%32;%mean(redDistEucl(:));
greenThreshold = 33;%mean(greenDistEucl(:));%40; %mean(greenDistEucl(:));
blueThreshold = 60;%mean(blueDistEucl(:));%25;%mean(blueDistEucl(:));

rgbDistEucl(:,:,1) = redDistEucl;
rgbDistEucl(:,:,2) = greenDistEucl;
rgbDistEucl(:,:,3) = blueDistEucl;
MaskRGB = rgbDistEucl(:,:,1)<redThreshold & rgbDistEucl(:,:,2) < greenThreshold & rgbDistEucl(:,:,3) < blueThreshold;

figure 
imshow(MaskRGB)


%% HSV

% % user determined mask
% I = rgb2hsv(imgFront);
% 
% % Define thresholds for channel 1 based on histogram settings
% channel1Min = 0.778;
% channel1Max = 0.046;
% 
% % Define thresholds for channel 2 based on histogram settings
% channel2Min = 0.106;
% channel2Max = 0.873;
% 
% % Define thresholds for channel 3 based on histogram settings
% channel3Min = 0.082;
% channel3Max = 0.671;
% 
% % Create mask based on chosen histogram thresholds
% BW = ( (I(:,:,1) >= channel1Min) | (I(:,:,1) <= channel1Max) ) & ...
%     (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
%     (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);


% convert to HSV
imgBottomHSV = double(rgb2hsv(imgBottom));
imgFrontHSV = double(rgb2hsv(imgFront));
%separate channels 
hChannelBottom = imgBottomHSV(:,:,1);
sChannelBottom = imgBottomHSV(:,:,2);
vChannelBottom = imgBottomHSV(:,:,3);

hChannelFront = imgFrontHSV(:,:,1);
sChannelFront = imgFrontHSV(:,:,2);
vChannelFront = imgFrontHSV(:,:,3);

% determine mean
hMeanBottom = mean(hChannelBottom(:));
sMeanBottom = mean(sChannelBottom(:));
vMeanBottom = mean(vChannelBottom(:));

hB = hMeanBottom*ones(size(hChannelFront,1),size(hChannelFront,2));
sB = sMeanBottom*ones(size(sChannelFront,1),size(sChannelFront,2));
vB = vMeanBottom*ones(size(vChannelFront,1),size(vChannelFront,2));

hDistEucl = sqrt((hChannelFront-hB).^2);hmin = min(hDistEucl(:));hmax = max(hDistEucl(:));
sDistEucl = sqrt((sChannelFront-sB).^2);smin = min(sDistEucl(:));smax = max(sDistEucl(:));
vDistEucl = sqrt((vChannelFront-vB).^2);vmin = min(vDistEucl(:));vmax = max(vDistEucl(:));


hsvDistEucl(:,:,1) = hDistEucl;
hsvDistEucl(:,:,2) = sDistEucl;
hsvDistEucl(:,:,3) = vDistEucl;
% set thresholds
hThreshold = mean(hDistEucl(:));%0.027;%0.03;%min(prctile(hDistEucl,95));%0.020;
sThreshold = mean(sDistEucl(:));%0.25;%0.36; %min(prctile(sDistEucl,95));%0.1; 
vThreshold = mean(vDistEucl(:));%0.45;%0.11; %min(prctile(vDistEucl,70));%0.1;
MaskHSV= hsvDistEucl(:,:,1)<=hThreshold & hsvDistEucl(:,:,2)<=sThreshold & hsvDistEucl(:,:,3)<=vThreshold;
figure
imshow(MaskHSV)

%% LAB
% convert to HSV
imgBottomLab = double(rgb2lab(imgBottom));
imgFrontLab = double(rgb2lab(imgFront));
%separate channels 
LChannelBottom = imgBottomLab(:,:,1);
aChannelBottom = imgBottomLab(:,:,2);
bChannelBottom = imgBottomLab(:,:,3);

LChannelFront = imgFrontLab(:,:,1);
aChannelFront = imgFrontLab(:,:,2);
bChannelFront = imgFrontLab(:,:,3);

% determine mean
LMeanBottom = mean(LChannelBottom(:));
aMeanBottom = mean(aChannelBottom(:));
bMeanBottom = mean(bChannelBottom(:));

LB = LMeanBottom*ones(size(LChannelFront,1),size(LChannelFront,2));
aB = aMeanBottom*ones(size(aChannelFront,1),size(aChannelFront,2));
bB = bMeanBottom*ones(size(bChannelFront,1),size(bChannelFront,2));

LDistEucl = sqrt((LChannelFront-LB).^2);Lmin = min(LDistEucl(:));Lmax = max(LDistEucl(:));
aDistEucl = sqrt((aChannelFront-aB).^2);amin = min(aDistEucl(:));amax = max(aDistEucl(:));
bDistEucl = sqrt((bChannelFront-bB).^2);bmin = min(bDistEucl(:));bmax = max(bDistEucl(:));


LabDistEucl(:,:,1) = LDistEucl;
LabDistEucl(:,:,2) = aDistEucl;
LabDistEucl(:,:,3) = bDistEucl;
% set thresholds
LThreshold = mean(LDistEucl(:));%12;%0.03;%min(prctile(hDistEucl,95));%0.020;
aThreshold = mean(aDistEucl(:));%4%12;%0.36; %min(prctile(sDistEucl,95));%0.1; 
bThreshold =mean(bDistEucl(:));%11%8;%0.11; %min(prctile(vDistEucl,70));%0.1;
MaskLab= LabDistEucl(:,:,1)<=LThreshold & LabDistEucl(:,:,2)<=aThreshold & LabDistEucl(:,:,3)<=bThreshold;
figure
imshow(MaskLab)