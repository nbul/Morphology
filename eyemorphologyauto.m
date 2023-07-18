clear variables;
close all;

%% Open the file
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
files = dir(strcat(filedir,'/*', '.tif'));
cd(filedir);

%Folder to save information about cells
if exist([filedir, '/ROI'],'dir') == 0
    mkdir(filedir,'/ROI');
end
roi_dir = [filedir, '/ROI'];

%Folder to save information about cells
if exist([filedir, '/threshold'],'dir') == 0
    mkdir(filedir,'/threshold');
end
thr_dir = [filedir, '/threshold'];

roi_all = zeros(numel(files),4);

cutoff = [0.8,0.25];

features = zeros(numel(files),8);
Distancesall = zeros(1,2);
Ommatidia = zeros(1,3);
for loop = 1:numel(files)
    cd(filedir);

    Im = imread(files(loop).name);
    ImRes = imresize(Im,0.25);


    se = strel('diamond',5);
    ImTH = imtophat(ImRes,se);
    ImTH = rgb2gray(ImTH);
    ImTH = imbinarize(ImTH);
    [ro, co] = find(ImTH == 1);
    rowcolCoordinates = [mean(ro), mean(co)];
    cc=bwconncomp(ImTH);
    Centers=regionprops(cc,'Centroid');
    C = [Centers.Centroid];
    CX = C(1:2:end);
    CY = C(2:2:end);
    C2 = [CX', CY'];
    Distances2 = sqrt((C2(:,1)-rowcolCoordinates(2)).^2 + (C2(:,2)-rowcolCoordinates(1)).^2 );
    %histogram(Distances2, 30)
    

    %% Create ROI mask
    [nx,ny] = size(ImRes(:,:,1));
    [X,Y] = meshgrid(1:ny,1:nx);
    th = linspace(0,2*pi);
    pLogic = (Distances2 < quantile(Distances2, 0.6));
    
    Coordinates2 = C2(pLogic == 1,:);
    rowcolCoordinates2 = mean(C2(pLogic == 1,:));
    Distances3 = sqrt((Coordinates2(:,1)-rowcolCoordinates2(1)).^2 + (Coordinates2(:,2)-rowcolCoordinates2(2)).^2 );
    pLogic2 = (Distances3 < quantile(Distances3, 0.1));
    Radius = quantile(Distances3, 0.1);

    xc = rowcolCoordinates2(1)+Radius*cos(th);
    yc = rowcolCoordinates2(2)+Radius*sin(th);
    idx = inpolygon(X(:),Y(:),xc',yc);
    image1=figure;
    imshow(ImRes)
    hold on
    plot(xc,yc,'b','LineWidth',2)
    cd(roi_dir);
    image_filename = append(regexp(files(loop).name,'\d*','Match'),'_ROI.tif');
    print(image1, '-dtiff', '-r150', char(image_filename));
    close all

    ImRes3 = ImRes;
    for i= 1:1:3
        ImRes2 = ImRes(:,:,i);
        ImRes2(~idx) = 0;
        ImRes3(:,:,i) = ImRes2;
    end


    %% Crop RGB image
    ImRes4 = ImRes3;
    ImRes4( ~any(ImRes3(:,:,1),2), : ,:) = [];  %rows
    ImRes4( :, ~any(ImRes3(:,:,1),1),: ) = [];


    %% Threshold and ommatidia centres
    Thresh = adaptthresh(imadjust(rgb2gray(ImRes)), 'NeighborhoodSize', 5);
    Im_thresh = imbinarize(imadjust(rgb2gray(ImRes)), 0.03 + Thresh);
    Im_thresh2 = Im_thresh;
    Im_thresh2(~idx) = 0;
    Im_thresh2( ~any(ImRes3(:,:),2), :) = [];  %rows
    Im_thresh2( :, ~any(ImRes3(:,:,1),1),: ) = [];
    Im_thresh2 = bwareaopen(Im_thresh2,3);
    % se = strel('disk',1);
    % Im_thresh2 = imopen(Im_thresh2, se);

    cd(thr_dir);
    image_filename = append(regexp(files(loop).name,'\d*','Match'),'_mask.tif');
    imwrite(Im_thresh2, char(image_filename));

    comp = bwconncomp(Im_thresh2);
    Dots = regionprops(comp, 'Centroid');

    %% Distances and regularity
    Dots2 = cat(1, Dots.Centroid);
    DT = delaunayTriangulation(Dots2);
    EdgesDots = edges(DT);
    Distances = zeros(size(EdgesDots,1),1);
    for i=1:1:size(EdgesDots,1)
        Distances(i) = sqrt((Dots2(EdgesDots(i,1),1)-Dots2(EdgesDots(i,2),1))^2 +...
            (Dots2(EdgesDots(i,1),2)-Dots2(EdgesDots(i,2),2))^2);
    end
    features(loop,1:3) = [str2double(regexp(files(loop).name,'\d*','Match')); mean(Distances); std(Distances)];
    Distancesall = [Distancesall; loop*ones(size(Distances,1),1), Distances];
    %% Voronoi
    [V,r] = voronoiDiagram(DT);
    false = zeros(size(V,1),1);
    false(V(:,1)<0) = 1;
    false(V(:,2)<0) = 1;
    false(V(:,1)>size(Im_thresh2,1)) = 1;
    false(V(:,2)>size(Im_thresh2,2)) = 1;

    clear area vertices
    idx2 = find(false ==1);
    counter = 0;
    for i=1:numel(r)
        if isempty(intersect(idx2,r{i})) == 1
            coordinates = V(r{i},:);
            counter = counter+1;
            area(counter) = polyarea(coordinates(:,1),coordinates(:,2));
            vertices(counter) = numel(r{i});
        end
    end
    features(loop,4:8) = [mean(area), std(area), mean(vertices), std(vertices), median(vertices)];
    Ommatidia = [Ommatidia; loop*ones(size(area',1),1), area', vertices'];
end


cd(roi_dir);
roi2 = array2table(roi_all);
roi2.Properties.VariableNames = {'eye', 'X', 'Y','Radius'};
writetable(roi2,'ROI.xlsx','Sheet','ROI','WriteMode','overwritesheet');

slashIdx = strfind(filedir, '/');
cd(filedir(1:slashIdx(end-2)))
features = rmoutliers(features, 'mean');
features2 = array2table(features);
features2.Properties.VariableNames = {'eye', 'Distance', 'SD Distance',...
    'Area','SD Area', 'Neighbour', 'SD Neighbour', 'Median Neighbour'};

pathSegment = [filedir(slashIdx(end-2)+1:slashIdx(end-1)-1),'_',...
    filedir(slashIdx(end-1)+1:slashIdx(end)-1), '_',filedir(slashIdx(end)+1:end)];
writetable(features2,'Summary.xlsx','Sheet',pathSegment,'WriteMode','overwritesheet');

Distancesall(1,:) = [];
Distancesall = rmoutliers(Distancesall, 'mean');
Distancesall2 = array2table(Distancesall);
Distancesall2.Properties.VariableNames = {'eye', 'Distance'};
writetable(Distancesall2,'SummaryDistance.xlsx','Sheet',pathSegment,'WriteMode','overwritesheet');

Ommatidia(1,:) = [];
Ommatidia = rmoutliers(Ommatidia, 'mean');
Ommatidia2 = array2table(Ommatidia);
Ommatidia2.Properties.VariableNames = {'eye', 'Area','Edges'};
writetable(Ommatidia2,'SummaryOmmatidia.xlsx','Sheet',pathSegment,'WriteMode','overwritesheet');

cd(currdir);
% %% if ey-GAL4 main folder
% cd(filedir(1:slashIdx(end-1)))
% features = rmoutliers(features, 'mean');
% features2 = array2table(features);
% features2.Properties.VariableNames = {'eye', 'Distance', 'SD Distance',...
%     'Area','SD Area', 'Neighbour', 'SD Neighbour', 'Median Neighbour'};
%
% pathSegment = [filedir(slashIdx(end-1)+1:slashIdx(end)-1), '_',filedir(slashIdx(end)+1:end)];
% writetable(features2,'Summary.xlsx','Sheet',pathSegment,'WriteMode','overwritesheet');