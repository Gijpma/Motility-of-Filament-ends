function lengthout=GetLength(sizebw,IDfilament,thickness)
%% calculates the length of the filament from ID data of bwconncomp. first 
% this ID data is converted to an bw bounding box image. then it calculates
% the average filament thickness by thinning the BW image and
% extrapolating how many thining operations would be needed to eliminate
% all filaments. from this thickness and the total ID count the filament
% length is calculated.

[x,y]=ind2sub(sizebw,IDfilament);
x=x-min(x)+1;
y=y-min(y)+1;
image=zeros(length(x),length(y));
image(sub2ind(size(image),x,y))=1;
lengthout=length(IDfilament)/thickness;
