function [frames, FilEnd, SizeFrame, num, Threshold,A]=VideoProcess(dirout, file, param)
%% this function reads the full video and does initial processing to find 
% filament ends and store the associate filament end images.
tic
VideoObject = VideoReader([dirout file]);
vidFrames = readFrame(VideoObject);
toc     % timer to indicate time taken for first video load. long time means
        % video stored on slow device or memory issues
MaxFrames=VideoObject.Duration*30;

% all results written to results folder.
if ~exist([dirout '\results'])
    mkdir([dirout '\results']);
end

%detect any dark edges and remove (some recording devices add black bands
%to compensate for a non-standard video sensor)
if find(vidFrames(:,0.5*end,1),1,'last')== find(vidFrames(:,0.5*end+10,1),1,'last')
edgeX=[find(vidFrames(:,0.5*end,1),1) find(vidFrames(:,0.5*end,1),1,'last')];
edgeY=[find(vidFrames(0.5*end,:,1),1) find(vidFrames(0.5*end,:,1),1,'last')];
else
    edgeX=[1 length(vidFrames(:,1,1))];
    edgeY=[1 length(vidFrames(1,:,1))];
end

%Store resulting image size
SizeFrame(1)=(range(edgeX)+1);
SizeFrame(2)=(range(edgeY)+1);

%Store multiple images for frame averaging
OldFrameStore(:,:,1)=vidFrames(edgeX(1):edgeX(2),edgeY(1):edgeY(2),1);
for i=2:param.FrameAverage-1
    vidFrames = readFrame(VideoObject);
    OldFrameStore(:,:,i)=vidFrames(edgeX(1):edgeX(2),edgeY(1):edgeY(2),1);
end


% initialize black and white image
bw=zeros(size(vidFrames(edgeX(1):edgeX(2),edgeY(1):edgeY(2),1)));

%% Main file reading and processing loop
for frames=param.FrameAverage:MaxFrames     % Starting at frame param.FrameAverage to allow averaging of frames
    
    if mod(frames,30)==0        % Indicates progress
        toc
        disp([num2str(frames) ' of ' num2str(MaxFrames) ' file ' file])
        tic
    end
    
    % load video frames.
    vidFrames = readFrame(VideoObject);
    % Frame load and cleaning of image
    OldFrameStore(:,:,param.FrameAverage)=vidFrames(edgeX(1):edgeX(2),edgeY(1):edgeY(2),1);
    A=mean(OldFrameStore,3);  % param.FrameAverage frame average for noise reduction
    OldFrameStore(:,:,1:param.FrameAverage-1)=OldFrameStore(:,:,2:param.FrameAverage);
    A=imgaussfilt(double(A(:,:,1)),1);   % 2d Gaussian filter of image
    A=A-0.8*imgaussfilt(A,20);           % Remove uneven lighting by subtracting heavily blurred image
    A=A*255/max(max(A));                 % adjust contrast to use the full 8bit gray scale range (0 255)
    if frames==param.FrameAverage
        %% Find treshold optimum by getting total number of elements for different threshold values
        for i=1:250
            bwcompfind=bwconncomp(A>i);    %use bwconncomp, much faster than bwlabel and find
            num(i)=bwcompfind.NumObjects;
        end
            [~,b1]=max(num);                % Identify peak element number location
                                            % optimal threshold a fixed
                                            % fraction above the point at
                                            % which 300 elements are found.
            Threshold=(b1(1)+find(num(b1(1):end)<300,1))*(1+param.BrightAdjust);
            brightness=sum(sum(A));         % average brightness of each frame is calculated 
                                            % and the threshold is adjusted
                                            % relative to the brightness
                                            % here
            if max(num) < 300               % specific for no noise videos
                Threshold = param.DefaultThreshold;
            end
    end
    % remove stray pixels from noise, both isolated (bwareaopen) and
    % attached to filaments (bwmorph majority)
    bw=bwareaopen(bwmorph(double(A>Threshold*(sum(sum(A)))/brightness),'majority',10),param.MinSize);
    % For filament length calculation. bwlabelout is used to identify
    % whether two filament ends occur in same connected component, while
    % bwlabelcc is used to calculate the length of the filament (faster).
    % Thickness is calculated by thinning all connected components by 1
    % pixels, and subsequently calculating what share of all pixels was
    % removed this way. This indicates the average thickness of each
    % filament. 
    bwlabelout=bwlabel(bw); 
    bwlabelcc=bwconncomp(bw);
    thickness=2*(nnz(bw)/(nnz(bw)-nnz(bwmorph(bw,'thin',1))));
    % first bw image is saved to confirm proper thresholding.
    if frames==param.FrameAverage
        imwrite(bw,[dirout '\results\' file(1:end-4) '_bw.png'])
    end
    
    % Morphological thining of total image followed by finding of all
    % endpoints
    Skeleton=bwmorph(bw,'thin','inf');
    ends=bwmorph(Skeleton,'endpoints');
    branch=bwmorph(Skeleton,'branchpoints');
    [xb,yb]=find(branch);
    %store all ends
    [x,y]=find(ends);
    rawx=x;rawy=y;   %used for finding filament length only
    % Remove ends that are too close to eachother or a branchpoint for accurate resolution
    for i=1:length(x)
        neighbours=find(abs(x(i)-x)<param.MaxDistAllow & abs(y(i)-y)<param.MaxDistAllow);    %finds if two ends are too close together
        neighbours2branch=find(abs(x(i)-xb)<param.MaxDistAllow*2 & abs(y(i)-yb)<param.MaxDistAllow*2);    %finds if two ends are to close together
        if length(neighbours)==2 
            if abs(diff(bwlabelout(sub2ind(size(bwlabelout),x(neighbours),y(neighbours)))))>0
                ends(sub2ind(size(bwlabelout),x(neighbours),y(neighbours)))=0;
            else
                ends(sub2ind(size(bwlabelout),x(neighbours),y(neighbours)))=0;
                [xtemp,ytemp]=find(bwlabelout==bwlabelout(x(neighbours(1)),y(neighbours(1))));
                ends(round(mean(xtemp)),round(mean(ytemp)))=1;
            end
        elseif length(neighbours)>2
            ends(sub2ind(size(bwlabelout),x(neighbours),y(neighbours)))=0;
        end
        if neighbours2branch
            ends(sub2ind(size(bwlabelout),x(i),y(i)))=0;
        end
    end
    
    
    
    % Remove any ends that are too close to the edge of the frame
    ends(1:SizeFrame(1),1:1*param.Square)=0; ends(1:SizeFrame(1),SizeFrame(2)-1*param.Square:SizeFrame(2))=0; ends(1:1*param.Square,1:SizeFrame(2))=0; ends(SizeFrame(1)-1*param.Square:SizeFrame(1),1:SizeFrame(2))=0;
    [x,y]=find(ends);
    
    labelno=bwlabelout(sub2ind(size(bw),rawx,rawy));  % for finding filament length
    
    
    FrameOut(frames).x=x;   %Stored to associate filament ends with eachother between frames
    FrameOut(frames).y=y;
    
    % for frame param.FrameAverage initial values are set, all other frames add to this
    if frames>param.FrameAverage
        
        endsout=[];
        for i=1:length(FrameOut(frames).x)
            % find correlation between tips and tails of current frame with
            % previous frame. based on absolute distance. faster than using
            % built in "dist"
            tempdist(:)=1000;  % setting to a constant is much faster then clearing the variable
            for j=1:length(FrameOut(frames-1).x)
                tempdist(j)=sqrt((FrameOut(frames).x(i)...
                    -FrameOut(frames-1).x(j))^2+(FrameOut(frames).y(i)...
                    -FrameOut(frames-1).y(j))^2);
            end
            [minval,minloc]=min(tempdist);   % finds minimum distance
            
            %define tipimage as a square of param.Square width and height
            tipimage=bw(x(i)-0.5*param.Square:x(i)+0.5*param.Square,y(i)-0.5*param.Square:y(i)+0.5*param.Square);
            
            % detection of proximity to other filament through
            % identification of more than one high-pixel area. 
            PositiveAreas=bwlabel(tipimage);                                
            % labelling the edge of the square tip area to detect direction
            % of movement. 
            edgelabel=(frames+0.1)*double((bwareaopen(tipimage,5)));          % label all pixels as frames + 0.1 
            edgelabel(2:end-1,2:end-1)=edgelabel(2:end-1,2:end-1)*frames/(frames+0.1);  %label all interior pixels as frames 
            [xout,yout]=find(edgelabel>0);                                    %find coordinates of all high pixels
            FilPoints=sub2ind(size(bw),xout+x(i)-0.5*param.Square-1,yout+y(i)-0.5*param.Square-1);   %adjust found points for full size image
            FilPoints(:,2)=edgelabel(find(edgelabel>0));                        % store value of pixels
            
            
            % if filament detected in proximity of tip, do not store 
            if max(max(PositiveAreas))==1
                 % only write new data to filament end if current location
                 % is less than Param.maxdistallow away, and filament end
                 % has been present less than param.MaxLength, and if a
                 % filament end has not yet been written to
                if minval<param.MaxDistAllow & frames-FilEnd(FilIDOld(minloc)).firstframe<param.MaxLength & ~(FilEnd(FilIDOld(minloc)).Coord(end,3)==frames)
                    
                    FilEnd(FilIDOld(minloc)).Coord(end+1,1:3)=[x(i) y(i) frames];                       %update position
                    FilIDNew(i)=FilIDOld(minloc);                                                         %update ID
                    FilEnd(FilIDOld(minloc)).TipImage(end+1:end+length(FilPoints(:,1)),:)=...
                        FilPoints;                                                                  % adds the new points to Tip Images
                    FilEnd(FilIDOld(minloc)).lastframe=frames;                                    % store last frame filament end was found
                    % get filament length only if two ends found in one
                    % connected component.
                    if length(find(labelno==bwlabelout(x(i),y(i))))==2
                         FilEnd(FilIDOld(minloc)).FilLength(end+1)=GetLength(size(bw),bwlabelcc.PixelIdxList{bwlabelout(x(i),y(i))},thickness);
                    else
                        FilEnd(FilIDOld(minloc)).FilLength(end+1)=0;
                    end
                else                                                        % a new filament end is stored
                    FilEnd(end+1).Coord(1,1:3)=[x(i) y(i) frames];        
                    FilIDNew(i)=length(FilEnd);                           
                    FilEnd(end).TipImage=double(FilPoints);
                    FilEnd(end).firstframe=frames;
                    FilEnd(end).reject=0;
                    if length(find(labelno==bwlabelout(x(i),y(i))))==2
                         FilEnd(end).FilLength(1)=GetLength(size(bw),bwlabelcc.PixelIdxList{bwlabelout(x(i),y(i))},thickness);
                    else
                        FilEnd(end).FilLength(1)=0;
                    end
                end
            else            % in this case no filament end is stored. dummy values placed
                endsout(end+1)=i;
                FilIDNew(i)=-1000;
                FrameOut(frames).x(i)=-1000;
                FrameOut(frames).y(i)=-1000;
            end
           
            
        end
  
        FilIDOld=FilIDNew;
        FilIDNew=[];
    else                                                %goes here only the first run, so all filament ends are new.
        for i=1:length(FrameOut(param.FrameAverage).x)
            FilEnd(i).Coord(1,1:3)=[FrameOut(param.FrameAverage).x(i) FrameOut(param.FrameAverage).y(i) frames];
            FilIDOld(i)=i;
            [xout,yout]=find(bw(x(i)-0.5*param.Square:x(i)+0.5*param.Square,y(i)-0.5*param.Square:y(i)+0.5*param.Square));          %find high points around ends
            FilPoints=sub2ind(size(bw),xout+x(i)-0.5*param.Square-1,yout+y(i)-0.5*param.Square-1);                    %adjust found points for full size image
            FilEnd(i).TipImage=double([FilPoints ones(length(FilPoints),1)*frames]);
            FilEnd(i).firstframe=frames;
            FilEnd(i).reject=0;
            if length(find(labelno==bwlabelout(x(i),y(i))))==2
                FilEnd(i).FilLength(1)=GetLength(size(bw),bwlabelcc.PixelIdxList{bwlabelout(x(i),y(i))},thickness);
            else
                FilEnd(i).FilLength(1)=0;
            end
        end
    end
    
end

%% removes traces that existed less than 10 frames
for i=1:length(FilEnd)
    if nnz(FilEnd(i).Coord(:,1))<param.MinPresence %| length(find(FilEnd(i).FinalImage-floor(FilEnd(i).FinalImage)))>0.2*length(find(FilEnd(i).FinalImage))
        K(i)=i;
    end
end
FilEnd(nonzeros(K))=[];

% save Raw data
save([dirout '\results\' file(1:end-4) '_ANL_Raw.mat'],'FilEnd','param','dirout','file','frames','VideoNow','num','Threshold','A','bw')