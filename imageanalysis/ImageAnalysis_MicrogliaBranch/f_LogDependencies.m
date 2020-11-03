function [FileList] = f_LogDependencies(Name, SavePath)
%Collect all upstream scripts and save them in a folder named Dependencies
%   Name: Function or script name as string
%   SavePath: 
%   Example: f_LogDependencies('Workflow_RosellaLC3_BPF_Summary_20170313.m', 'S:\HCS_Platform\Data\JavierJarazo\Rosella\LC3\Analysis_20170323_143443')
    
    SavePath = [SavePath filesep 'Dependencies'];
    mkdir(SavePath);
    
    %% Scan recursively %%
    FileList = matlab.codetools.requiredFilesAndProducts(Name)';
    %%%%%%%%%%%%%%%%%%%%%%
    FileCount = size(FileList,1);

    for i = 1:FileCount
        FileThis = regexp(FileList{i}, '\\.*\\(.*)', 'tokens'); FileThis = FileThis{:}{:};
        SavePathThisFile = [SavePath, filesep, FileThis];
        copyfile(FileList{i}, SavePathThisFile);
    end
    
end

