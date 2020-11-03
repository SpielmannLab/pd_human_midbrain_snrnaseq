
function [FieldMask] = f_Create_Mask_from_ObjectList_Pixel_IDX(ObjectList, IndexVarName, Image)
% Author: Paul Antony 20140717
% For an objectlist in table format, let IndexVarName be an object related pixel index.
% This function creates a mask which represents all objects of this objectlist
%   ObjectList: Name of the input table
%   IndexVarName: Variable name for index vectors provided as 'string'


%% Initialize mask
FieldMask = zeros(size(Image));
ListSize = size(ObjectList);

PixelData = eval(['ObjectList.', IndexVarName]);

    for i = 1:ListSize
        try
            PixelDataObject = PixelData{i};
            FieldMask(PixelDataObject) = 1;
        catch
            continue
        end
    end

end

