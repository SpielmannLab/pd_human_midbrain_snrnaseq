%% Analysis of neurons
%% Prepare working space
clear
clc
%% Document script
files = dirrec('Z:\groups\gruenewald\semra.smajic\20200924_SN\2020-09-24', '.czi');
files = files';
filesTable = cell2table(files);
MainPath = ['S:\HCS_Platform\Data\SemraSmajic\ForPaper\FinalSN\Analysis_C045_SN_', datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(MainPath)
PreviewPath = [MainPath, '\Previews'];
TablesPath = [MainPath, '\Tables'];
mkdir(PreviewPath)
mkdir(TablesPath)
FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[MainPath, '\', FileNameShort]);
FileNameAndLocation = mfilename('fullpath');
currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup); 
f_LogDependencies(FileNameShort, MainPath); 

%% Loading data

for f = 2 : size(files, 1)

    
    BfCube = bfopen(files{f});
    
    omeMeta = BfCube{1, 4};
    InfoTableThis = table();
    AreaName = regexp(files{f}, '\\.*\\(.*).czi', 'tokens');
    InfoTableThis.AreaName= string(AreaName{:}{:});   
    InfoTableThis.ImSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    InfoTableThis.ImSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    InfoTableThis.PixellSizeX = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in µm
    InfoTableThis.PixelSizeY = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in µm
    InfoTableThis.Totaltiles = size(BfCube,1);
    ThisPreviewPath = [PreviewPath, '\', AreaName{:}{:}];
    mkdir(ThisPreviewPath)
    
    parfor t= 1:size(BfCube, 1)
    ch1 = BfCube{t,1}{1,1}; % Hoechst
    ch2 = BfCube{t,1}{2,1}; % 488-CADPS2
    ch3 = BfCube{t,1}{3,1}; % 568-MAP2
    ch4 = BfCube{t,1}{4,1}; % 647-TH
    %ch5 = BfCube{t,1}{5,1}; % Brightfield_1
    ch6 = BfCube{t,1}{5,1}; % Brightfield_2  %
    it(ch4)
    
    ObjectsAllTiles{f,t} = ImageAnalysis_Neuron_seg(ch1,ch2,ch3,ch4,ch6,InfoTableThis, ThisPreviewPath,t,AreaName);
    
    
    end

end 
ObjectsAllTiles = vertcat(ObjectsAllTiles{:});

% Writa a table with all samples in a loop
writetable(ObjectsAllTiles, [TablesPath, filesep, 'All_sections.csv'], 'WriteVariableNames', true); % Saving as comma separated file
