clc
clear
% Specify the folder where the files live.
inputFolder = '...';
outputFolder = ['...', datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(outputFolder)

PreviewPath = [outputFolder, '\Previews'];
TablesPath = [outputFolder, '\Tables'];
mkdir(PreviewPath)
mkdir(TablesPath)

FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[outputFolder, '\', FileNameShort]);
FileNameAndLocation = mfilename('fullpath');
currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup); 
f_LogDependencies(FileNameShort, outputFolder);

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(inputFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nWill you please specify a new folder.', inputFolder);
    uiwait(warndlg(errorMessage));
    inputFolder = uigetdir(); % Ask for a new one.
    if inputFolder == 0
         % cancelled by user
         return;
    end
end
% Get a list of all files in the folder ('*.png') or subfolders ('**/*.png') with the desired file name pattern.
filePattern = fullfile(inputFolder, '**/*AF488.tif'); % Change to the pattern you need.
myFiles = dir(filePattern);

for k = 1:numel(myFiles)
  FileName = myFiles(k).name;
  sampleName = erase(FileName, '-TH647_PLP546_GFAP488_DAPI405'); 
  sampleName = erase(sampleName, '_AF488.tif');
  fullFileName = fullfile(myFiles(k).folder, FileName);
  rawImage = imread(fullFileName);   % it(rawImage)
  ThisPreviewPath = [PreviewPath, '\', sampleName];
  mkdir(ThisPreviewPath)

  InfoTableThis = table();
  InfoTableThis.sampleName= string(sampleName);

%% Image Analysis
  grayImage = rgb2gray(rawImage);   % it(grayImage)
 
  CellMask = grayImage > 10;
  CellMask = bwareaopen(CellMask, 200); 
  [CellMask, CellCount] = bwlabeln(CellMask);      
  CellMask = imfill(CellMask,'holes');
  
  ImObjects = regionprops('table', CellMask, grayImage, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
  fullOutputFileName = fullfile(myFiles(k).name);
  

  %%  Skeletonize microglial body
   BodyMask = imopen(CellMask, strel('disk', 5));  % it(BodyMask)   it(Branchpoints)  
   Skeleton = Skeleton3D(CellMask); %imtool(CellMask + Skeleton, [])
   Skeleton = Skeleton & ~BodyMask;  %imtool(Skeleton)
   Branchpoints = bwmorph(Skeleton, 'branchpoints'); %imtool(BodyMask + Skeleton + Branchpoints) 
   Endpoints = bwmorph(Skeleton, 'endpoints'); %imtool(CellMask + Skeleton + Branchpoints + (2*Endpoints), [])
   C = imoverlay(BodyMask, Branchpoints,[1 0 1]);
   %imtool(C)
  
   %% 

   for i = 1:CellCount
       imThisCell = CellMask == i;  % it(imThisCell)
       SkeletonThisCell = Skeleton .* imThisCell; % 
       SkelByArea = sum(SkeletonThisCell(:)) / sum(imThisCell(:));
       ImObjects(i, 'SkeletonArea') = {sum(SkeletonThisCell(:))};
       ImObjects(i, 'SkelByArea') = {SkelByArea};
   end
        
  %% Making final table
    Objectscombined =ImObjects(:, [1:6, 8, 10:14]);
  
    SavePathForObjects = [TablesPath, filesep, sampleName, '_', 'objects.csv'];
    writetable(Objectscombined, SavePathForObjects); % Saving as comma separated file

          
 
%% previews
  ImMaskPerimeter = bwperim(CellMask);
  ImMaskPerimeter = imdilate(ImMaskPerimeter, true(0.5));
  overlay = imoverlay(rawImage, ImMaskPerimeter, [1 1 1]); % imshow(overlay)
  
  % Save preview
  SavePathMaskPreview = [PreviewPath, filesep, sampleName, '_', '_IBA1_mask.png'];
  imwrite(overlay, SavePathMaskPreview)
  
  % Save raw image
  SavePathRawPreview = [PreviewPath, filesep, sampleName, '_', '_IBA1_raw.png'];
  imwrite(rawImage, SavePathRawPreview)

  % Save branched image
  SavePathRawPreview = [PreviewPath, filesep, sampleName, '_', '_IBA1_branch.png'];
  imwrite(C, SavePathRawPreview)
  
 
 end

