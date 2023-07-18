clear variables;
close all;

%% Open the file
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
files = dir(strcat(filedir,'/*', '.tif'));
cd(filedir);

cutoff = [0.8,0.25];

for loop = 1:numel(files)

    Im = imread([num2str(loop), '.tif']);
    ImRes = imresize(Im,0.25);
    se = strel('diamond',5);
    ImTH = imtophat(ImRes,se);
    ImTH = rgb2gray(ImTH);
    ImTH = imbinarize(ImTH);
    [r, c] = find(ImTH == 1);
    rowcolCoordinates = [mean(r), mean(c)];
    cc=bwconncomp(ImTH);
    Centers=regionprops(cc,'Centroid');
    C = [Centers.Centroid];
    CX = C(1:2:end);
    CY = C(2:2:end);
    C2 = [CX', CY'];
    Distances = sqrt((C2(:,1)-rowcolCoordinates(2)).^2 + (C2(:,2)-rowcolCoordinates(1)).^2 );
    pLogic = (Distances < quantile(Distances, cutoff(1)));
    % histogram(Distances, 30)
    % hold on;
    % line([quantile(Distances, 0.8), quantile(Distances, 0.8)], ylim, 'LineWidth', 2, 'Color', 'r');
    imshow(ones(size(ImTH)))
    hold on
    scatter(C2(pLogic == 0,1),C2(pLogic == 0, 2),8,'b','filled')
    hold on;
    scatter(C2(pLogic == 1,1),C2(pLogic == 1, 2),8,'r','filled')
    pLogic2 = (Distances < quantile(Distances, cutoff(2)));
    imshow(ImRes)
    R = quantile(Distances, cutoff(2));
    th = linspace(0,2*pi) ;
    xc = rowcolCoordinates(2)+R*cos(th) ;
    yc = rowcolCoordinates(1)+R*sin(th) ;
    hold on
    plot(xc,yc,'b','LineWidth',2)
    image_filename = append(regexp(files(loop).name,'\d*','Match'),'_ROI.tif');
    print(image1, '-dtiff', '-r150', char(image_filename));

    [nx,ny] = size(ImTH);
    [X,Y] = meshgrid(1:ny,1:nx);
    %imshow(ImRes) ;
    %hold on
    R = quantile(Distances, cutoff(2));
    th = linspace(0,2*pi);
    xc = rowcolCoordinates(2)+R*cos(th);
    yc = rowcolCoordinates(1)+R*sin(th);
    idx = inpolygon(X(:),Y(:),xc',yc);
    
    ImRes3 = ImRes;
    for i= 1:1:3
        ImRes2 = ImRes(:,:,i);
        ImRes2(~idx) = 0;
        ImRes3(:,:,i) = ImRes2;
    end
   
    ImRes4 = ImRes3;
    ImRes4( ~any(ImRes3(:,:,1),2), : ,:) = [];  %rows
    ImRes4( :, ~any(ImRes3(:,:,1),1),: ) = [];
    
    
    Thresh = adaptthresh(imadjust(rgb2gray(ImRes)), 'NeighborhoodSize', 5);
    Im_thresh = imbinarize(imadjust(rgb2gray(ImRes)), 0.002 + Thresh); 
    Im_thresh2 = Im_thresh;
    Im_thresh2(~idx) = 0;
    Im_thresh2( ~any(ImRes3(:,:),2), :) = [];  %rows
    Im_thresh2( :, ~any(ImRes3(:,:,1),1),: ) = [];
    
    se = strel('disk',1);
    Im_thresh2 = imopen(Im_thresh2, se);
    nmask = fillHull(nmask); #display(nmask)
    nmask = bwlabel(nmask)

end