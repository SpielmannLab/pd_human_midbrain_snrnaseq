function [BarMask, BarCenter] = f_barMask(length_um, umPerPixel, imSize, position_x, position_y, height)
%Create a mask corresponding to a scalebar:
    %   length_um defines the length of the scalebar in um
    %   umPerPixel corresponds to the pixelsize in um per pixel
    %   imSize corresponds to the size of the image. Use size(ImageOfInterest) to assign this value
    %   position corresponds to the bottom left pixel coordinate of the scalebar
    %       x,
    %       y,
    %   height corresponds to the hight of the scalebar

Image = zeros(imSize(1), imSize(2));
BarPixelLength = round(length_um / umPerPixel);
%Image((position_x : position_x + BarPixelLength), (position_y : position_y + hight)) = 1;
Image((position_x : position_x + height), (position_y : position_y + BarPixelLength)) = 1;
BarMask = Image;
BarMask = (BarMask ~= 0);
BarCenterX = (position_x + (position_x + height)) / 2;
BarCenterY = (position_y + (position_y + BarPixelLength)) /2;
BarCenter = [BarCenterX, BarCenterY];
 

end

