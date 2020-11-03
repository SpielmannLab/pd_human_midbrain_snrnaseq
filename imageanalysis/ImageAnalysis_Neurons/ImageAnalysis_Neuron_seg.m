function [Objectscombined] = ImageAnalysis_Neuron_seg(ch1,ch2,ch3,ch4,ch6,InfoTableThis, ThisPreviewPath,t,AreaName)  
	
    %% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',1, 1)); 
    NucleiMask = NucleiBlurred > 300;
    NucleiLN = bwlabeln(NucleiMask);
    NucleiObjects = regionprops('table', NucleiLN, ch1, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    
    NucleiObjects = NucleiObjects(NucleiObjects.Area < 15000,:);
    NucleiObjects = NucleiObjects(NucleiObjects.Area > 50,:);
    NucleiObjects = NucleiObjects(NucleiObjects.Eccentricity < 0.85,:);
    NucleiMask = f_Create_Mask_from_ObjectList_Pixel_IDX(NucleiObjects,'PixelIdxList', ch1);
       
    NucleiObjects = NucleiObjects(:,[1:3,10:12]);
    NucleiObjects.ObjectType = repmat("Nuclei",height(NucleiObjects),1);
    NucleiObjects.Tilenumber = repmat(t,height(NucleiObjects),1);       
  
	
    %% Segment TH
    THBlurred = imfilter(ch4, fspecial('gaussian',2, 1));
    
    THmean = mean2(THBlurred);
    % Calculate meanInt of an image and keep whatever is 1.2 fold higher
    THMask = THBlurred > THmean * 1.2; 
    THMask = bwareaopen(THMask, 200);
    
    % Segment with watershed
    D_TH = -bwdist(~THMask);
    % Compute the watershed transform of D image
    % Segment the image by watershed ridge lines (Ld ==0) by setting their
    % pixel value to 0 -> apply it to bw image
    mask2_TH = imextendedmin(D_TH,2);   
    D2_TH = imimposemin(D_TH,mask2_TH);
    Ld2_TH = watershed(D2_TH);
    mask3_TH = THMask;
    mask3_TH(Ld2_TH == 0) = 0; 
    %% Soma TH Objects
    THLN = bwlabeln(mask3_TH);
    THObjects = regionprops('table', THLN, ch4, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    
    % 2 filters
    Filtered_1 = THObjects(THObjects.Area > 800 & THObjects.Area < 20000,:);
    Filtered_1 = Filtered_1(Filtered_1.MajorAxisLength > 70,:);
    Filtered_1 = Filtered_1(Filtered_1.MinorAxisLength > 50,:);
   
    TH_Mask_Filt_1 = f_Create_Mask_from_ObjectList_Pixel_IDX(Filtered_1,'PixelIdxList', ch4); 
    
    THSomaMask_final = TH_Mask_Filt_1; 
    THSomaObjects_final = vertcat(Filtered_1);
    THSomaObjects_final = THSomaObjects_final(:,[1:3,10:12]);
    THSomaObjects_final = unique(THSomaObjects_final,'rows');
    THSomaObjects_final.ObjectType = repmat("THpos",height(THSomaObjects_final),1);
    THSomaObjects_final.Tilenumber = repmat(t,height(THSomaObjects_final),1);
    
    
    %% Segment MAP2
    MAP2Blurred= imfilter(ch3, fspecial('gaussian',10, 1)); 
    
    % Make a mask to be removed
    RemMAP2_1 = MAP2Blurred < 300;    
    RemMAP2_2 = MAP2Blurred > 700;   
    RemMAP2_1 = imdilate(RemMAP2_1,strel('disk',3));  
    RemMAP2_2 = imdilate(RemMAP2_2,strel('disk',3));  
    
	MAP2Mask = MAP2Blurred > 300; 
    MAP2Mask = bwareaopen(MAP2Mask, 2000, 4); 
	
    MAP2SomaMask = MAP2Mask & ~RemMAP2_1;  
    MAP2SomaMask = MAP2SomaMask & ~RemMAP2_2;
    MAP2SomaMask = imfill(MAP2SomaMask, 'holes');
    MAP2SomaMask =  bwareaopen(MAP2SomaMask, 2000, 4);
    
    % Smaller objects are connected to the main object (cell soma)
    % Segment with watershed
    D = -bwdist(~MAP2SomaMask); 
    % Compute the watershed transform of D image
    mask2 = imextendedmin(D,2);    
    D2 = imimposemin(D,mask2);
    Ld2 = watershed(D2);
    mask3 = MAP2SomaMask;
    mask3(Ld2 == 0) = 0;     
    
    % After segmentation, we remove very small
    % objects which were attached to the main object by Area in regionprops
  
    MAP2LN = bwlabeln(mask3);
    
    MAP2Objects = regionprops('table', MAP2LN, ch3, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    MAP2FinalObjects = MAP2Objects(MAP2Objects.Eccentricity < 0.85,:);
    MAP2FinalObjects = MAP2FinalObjects(MAP2FinalObjects.Area > 1500,:);
    MAP2FinalObjects = MAP2FinalObjects(MAP2FinalObjects.Area < 20000,:);
    MAP2FinalObjects = MAP2FinalObjects(MAP2FinalObjects.MajorAxisLength > 30,:);
    MAP2FinalObjects = MAP2FinalObjects(MAP2FinalObjects.MinorAxisLength > 10,:);
    
    MAP2FinalSomaMask = f_Create_Mask_from_ObjectList_Pixel_IDX(MAP2FinalObjects,'PixelIdxList', ch3); 
    MAP2FinalObjects = MAP2FinalObjects(:,[1:3,10:12]);
    MAP2FinalObjects.ObjectType = repmat("MAP2pos",height(MAP2FinalObjects),1);
    MAP2FinalObjects.Tilenumber = repmat(t,height(MAP2FinalObjects),1);
    %% CADPS2 intensity in all neurons (MAP2+)  - raw ch2 intensity
    CADPS2MaskNEURO = bwlabeln(MAP2FinalSomaMask);  
    CADPS2ObjectsNEURO = regionprops('table', CADPS2MaskNEURO, ch2, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    
    CADPS2ObjectsNEURO = CADPS2ObjectsNEURO(:,[1:3,10:12]);
    CADPS2ObjectsNEURO.ObjectType = repmat("CADPS2inMAP2",height(CADPS2ObjectsNEURO),1);
    CADPS2ObjectsNEURO.Tilenumber = repmat(t,height(CADPS2ObjectsNEURO),1);
    %% CADPS2 intensity in DANs (TH+)  - raw ch2 intensity
    CADPS2MaskDAN = bwlabeln(THSomaMask_final);  
    CADPS2ObjectsDAN = regionprops('table', CADPS2MaskDAN, ch2, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
   
    CADPS2ObjectsDAN = CADPS2ObjectsDAN(:,[1:3,10:12]);
    CADPS2ObjectsDAN.ObjectType = repmat("CADPS2inTH",height(CADPS2ObjectsDAN),1);
    CADPS2ObjectsDAN.Tilenumber = repmat(t,height(CADPS2ObjectsDAN),1);
    %% Segment Brightfield
    BFBlurred= imfilter(ch6, fspecial('gaussian',10, 8)); 
    BFMask = BFBlurred < 1700;
    NeuromelMask =  bwareaopen(BFMask, 300); 

    NeuromelLN = bwlabeln(NeuromelMask);
    NeuromelObjects = regionprops('table', NeuromelLN, ch6, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    NMFinalObjects = NeuromelObjects(NeuromelObjects.Eccentricity < 0.85,:);
    NMFinalObjects = NMFinalObjects(NMFinalObjects.Area > 300 & NMFinalObjects.Area < 50000,:);

    NMFinalMask = f_Create_Mask_from_ObjectList_Pixel_IDX(NMFinalObjects,'PixelIdxList', ch6);  
    NMFinalObjects = NMFinalObjects(:,[1:3,10:12]);
    NMFinalObjects.ObjectType = repmat("Neuromel",height(NMFinalObjects),1);
    NMFinalObjects.Tilenumber = repmat(t,height(NMFinalObjects),1);       
    
    %% Expansion of Neuromel for colocalization
    NeuromelMaskdil = imdilate(NMFinalMask,strel('disk',10));
   
    %% MAP2 positive TH negative
    CatMAP2posTHneg = MAP2FinalSomaMask & ~THSomaMask_final; 
    CatMAP2posTHneg = bwareaopen(CatMAP2posTHneg, 250); 
    
    CatMAP2posTHnegdil = imdilate(CatMAP2posTHneg,strel('disk',10));
    
    MAP2posTHnegLN = bwlabeln(CatMAP2posTHneg);
    MAP2posTHnegObjects = regionprops('table', MAP2posTHnegLN, ch3, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    MAP2posTHnegObjects = MAP2posTHnegObjects(:,[1:3,10:12]);
    MAP2posTHnegObjects.ObjectType = repmat("MAP2posTHneg",height(MAP2posTHnegObjects),1);
    MAP2posTHnegObjects.Tilenumber = repmat(t,height(MAP2posTHnegObjects),1);
    
    %% MAP2 positive TH negative Neuromel pos
    CatMAP2posTHnegNeuromelpos = CatMAP2posTHnegdil & NeuromelMaskdil; %
	CatMAP2posTHnegNeuromelpos = bwareaopen(CatMAP2posTHnegNeuromelpos, 250);
    
    MAP2posTHnegNeuromelposLN = bwlabeln(CatMAP2posTHnegNeuromelpos); 
    MAP2posTHnegNeuromelposObjects = regionprops('table', MAP2posTHnegNeuromelposLN, ch3, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    MAP2posTHnegNeuromelposObjects = MAP2posTHnegNeuromelposObjects(:,[1:3,10:12]);
    MAP2posTHnegNeuromelposObjects.ObjectType = repmat("MAP2posTHnegNeuromelpos",height(MAP2posTHnegNeuromelposObjects),1);
    MAP2posTHnegNeuromelposObjects.Tilenumber = repmat(t,height(MAP2posTHnegNeuromelposObjects),1);
    %% CADPS2 intensity in  MAP2 positive, TH negative, Neuromel pos   - raw ch2 intensity
    CADPS2MaskMAP2NM = bwlabeln(MAP2posTHnegNeuromelposLN); 
    CADPS2MaskMAP2NM = regionprops('table', CADPS2MaskDAN, ch2, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
   
    CADPS2MaskMAP2NM = CADPS2MaskMAP2NM(:,[1:3,10:12]);
    CADPS2MaskMAP2NM.ObjectType = repmat("CADPS2inMapAndNM",height(CADPS2MaskMAP2NM),1);
    CADPS2MaskMAP2NM.Tilenumber = repmat(t,height(CADPS2MaskMAP2NM),1);
	%% CADPS2 intensity in MAP2 positive, TH positive, Neuromel pos   - raw ch2 intensity
    CatTHposNeuromelpos = MAP2FinalSomaMask & THSomaMask_final & NeuromelMaskdil; 
    THposNeuromelposLN = bwlabeln(CatTHposNeuromelpos);
    THposNeuromelposObjects = regionprops('table', THposNeuromelposLN, ch2, {'Area','ConvexArea','Solidity','PixelValues', 'PixelIdxList','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity'});
    THposNeuromelposObjects = THposNeuromelposObjects(:,[1:3,10:12]);
    THposNeuromelposObjects.ObjectType = repmat("CADPS2inTHposNeuromelpos",height(THposNeuromelposObjects),1);
    THposNeuromelposObjects.Tilenumber = repmat(t,height(THposNeuromelposObjects),1);
    
    %% Making final table
    Objectscombined =vertcat(NucleiObjects,THSomaObjects_final,CADPS2ObjectsDAN,MAP2FinalObjects,CADPS2ObjectsNEURO,NMFinalObjects,MAP2posTHnegObjects,MAP2posTHnegNeuromelposObjects,CADPS2MaskMAP2NM,THposNeuromelposObjects);
    Objectscombined = horzcat(repmat(InfoTableThis,height(Objectscombined),1),Objectscombined);  
              
    %% 2D previews
    
    % Scalebar
    imSize = size(ch3);
    PixelSizeX = InfoTableThis.ImSizeX;
	[BarMask, ~] = f_barMask(20, PixelSizeX, imSize, imSize(1)-50, 75, 10);

    RGB_1 = cat(3, imadjust(ch3,[0 0.01], [0 1]), imadjust(ch4), imadjust(ch1));
    RGB_1 = imoverlay(RGB_1, BarMask, [1 1 1]);

    RGB_2 = cat(3, imadjust(ch3,[0 0.01], [0 1]), imadjust(ch4), imadjust(ch6));
	RGB_2 = imoverlay(RGB_2, BarMask, [1 1 1]); % imtool(RGB_2)
  
    MAP2MaskPreview = imoverlay(imadjust(ch3, [0 0.01], [0 1]), bwperim(MAP2FinalSomaMask), [1 0 0]); %imtool(MAP2MaskPreview )
    MAP2RawPreview = imadjust(ch3, [0 0.01], [0 1]); %imtool(MAP2RawPreview)
	
	THMaskPreview = imoverlay(imadjust(ch4), bwperim(THSomaMask_final), [0 1 0]); %imtool(THMaskPreview )
    THRawPreview = imadjust(ch4);
   
	NucleiMaskPreview = imoverlay(imadjust(ch1), bwperim(NucleiMask), [1 0 0]); %imtool(NucleiMaskPreview )
    NucleiRawPreview = imadjust(ch1);
    
    BFMaskPreview = imoverlay(imadjust(ch6), bwperim(NMFinalMask), [1 0 0]); %imtool(BFMaskPreview )
    BFRawPreview = imadjust(ch6); %imtool(BFRawPreview )
    
    
    SavePathMAP2MaskPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_MAP2_mask.png'];
    SavePathMAP2RawPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_MAP2_raw.png'];
    SavePathTHMaskPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_TH_mask.png'];
    SavePathTHRawPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_TH_raw.png'];
	SavePathNucleiMaskPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_Nuclei_mask.png'];
    SavePathNucleiRawPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_Nuclei_raw.png'];
	SavePathBFMaskPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_BF_mask.png'];
    SavePathBFRawPreview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_BF_raw.png'];
    
	SavePathRGB_1Preview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_RGB_1.png'];
    SavePathRGB_2Preview = [ThisPreviewPath, filesep, AreaName{:}{:},'_', num2str(t), '_RGB_2.png'];
    
    imwrite(MAP2MaskPreview, SavePathMAP2MaskPreview)
    imwrite(MAP2RawPreview,  SavePathMAP2RawPreview)
	imwrite(NucleiMaskPreview, SavePathNucleiMaskPreview)
    imwrite(NucleiRawPreview,  SavePathNucleiRawPreview)
	imwrite(THMaskPreview, SavePathTHMaskPreview)
    imwrite(THRawPreview,  SavePathTHRawPreview)
    imwrite(BFMaskPreview, SavePathBFMaskPreview)
    imwrite(BFRawPreview,  SavePathBFRawPreview)
    imwrite(RGB_1, SavePathRGB_1Preview)
	imwrite(RGB_2, SavePathRGB_2Preview)
    
       
end