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
f_LogDependencies(FileNameShort, MainPath); 

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(inputFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', inputFolder);
    uiwait(warndlg(errorMessage));
    inputFolder = uigetdir(); % Ask for a new one.
    if inputFolder == 0
         
         return;
    end
end
% Get a list of all files in the folder ('*.png') or subfolders ('**/*.png') with the desired file name pattern.
filePattern = fullfile(inputFolder, '**/*AF488.tif'); % Change to the pattern you need, depending on which chanel you analyze
myFiles = dir(filePattern);

for k = 1:numel(myFiles)
  FileName = myFiles(k).name;
  sampleName = erase(FileName, '-TH647_PLP546_GFAP488_DAPI405'); 
  sampleName = erase(sampleName, '_AF488.tif');
  fullFileName = fullfile(myFiles(k).folder, FileName);
  rawImage = imread(fullFileName);   
  ThisPreviewPath = [PreviewPath, '\', sampleName];
  mkdir(ThisPreviewPath)

  InfoTableThis = table();
  InfoTableThis.sampleName= string(sampleName);

%% Image Analysis
  grayImage = rgb2gray(rawImage);   
  ImMask = grayImage > 10;    %IBA1>10; PLP1>5; GFAP>10;
  ImMask = bwlabeln(ImMask);     
  ImMask = imfill(ImMask,'holes');
     
  ImObjects = regionprops('table', ImMask, grayImage, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
  fullOutputFileName = fullfile(myFiles(k).name);
  
  ImObjects = ImObjects(ImObjects.Area > 200,:);
  ImMask = f_Create_Mask_from_ObjectList_Pixel_IDX(ImObjects,'PixelIdxList', grayImage);
      
  %% Making final table
    Objectscombined =ImObjects(:, [1:6, 8, 10:12]);
   
    SavePathForObjects = [TablesPath, filesep, sampleName, '_', 'objects.csv'];
    writetable(Objectscombined, SavePathForObjects); % Saving as comma separated file

   
 %% previews
   ImMaskPerimeter = bwperim(ImMask);
   ImMaskPerimeter = imdilate(ImMaskPerimeter, true(0.5));
   overlay = imoverlay(rawImage, ImMaskPerimeter, [1 1 1]); 
   
   % Save preview - rename as needed
   SavePathMaskPreview = [PreviewPath, filesep, sampleName, '_', '_IBA1_mask.png'];
   imwrite(overlay, SavePathMaskPreview)
   
  % Save raw image - rename as needed
   SavePathRawPreview = [PreviewPath, filesep, sampleName, '_', '_IBA1_raw.png'];
   imwrite(rawImage, SavePathRawPreview)

 
 end

