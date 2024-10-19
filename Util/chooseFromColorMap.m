function [Color, PreviouslyUsedColorIndexes] = chooseFromColorMap(PreviouslyUsedColorIndexes, NColor)

   ColorMap =[128 0 0
              230 25 75
              250 190 212
              170 110 40
              245 130 48
              255 125 180
              128 128 0
              255 225 25
%               255 250 200
              210 245 60 %10
              60 180 75
              170 255 195
              0 128 128 
              70 240 240
              0 0 128
              0 130 200
              145 30 180
              220 190 255
              240 50 230
%               0 0 0
              128 128 128]*1/255;
          
    Indexes = 1:length(ColorMap);
    Color = zeros(NColor, 3);
    ColorMap(Indexes(PreviouslyUsedColorIndexes),:)= [];
    for i = 1:NColor
        ColorIndex = ceil(rand(1)*(size(ColorMap, 1)));
        if ColorIndex == 0
            Color(i, :) = rand(1, 3);
        else
            Color(i, :) = ColorMap(ColorIndex, :);
            PreviouslyUsedColorIndexes = [PreviouslyUsedColorIndexes ColorIndex];
            ColorMap(Indexes(PreviouslyUsedColorIndexes(end)),:)= [];
        end

    end
end