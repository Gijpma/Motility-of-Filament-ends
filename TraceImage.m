function [FilEnd]=TraceImage(FilEnd, SizeFrame)
%% generates trace images

for i=1:length(FilEnd)
    [x,y]=ind2sub(SizeFrame, FilEnd(i).TipImage(:,1));            % first reduce image size to just the trace boundingbox
    FilEnd(i).historyshift=[min(y) min(x)];                       %and store where this boundingbox is
    x=x-min(x)+1;
    y=y-min(y)+1;
    indxy=sub2ind([max(x) max(y)],x,y);
    ImageTail=(zeros(max(x),max(y)));
    ImageTail(indxy)= FilEnd(i).TipImage(:,2);                % generate trace image with maximum frame values
    ImageTip=(zeros(max(x),max(y)));% store
    ImageTip(indxy(end:-1:1))= FilEnd(i).TipImage(end:-1:1,2);                 % generate image with minimum frame falues
    OvertakeTip=double(bwmorph(ImageTip,'thin',2)).*ImageTip;
    OvertakeTail=double(bwmorph(ImageTail,'thin',2)).*ImageTail;
    clear TailOut
    % detect Tail_Tip Overtake
    for j=FilEnd(i).firstframe:FilEnd(i).lastframe
        TailOut(j-FilEnd(i).firstframe+1,1)=j;
        TailOut(j-FilEnd(i).firstframe+1,2)=length(find(OvertakeTail==j+0.1));
        TailOut(j-FilEnd(i).firstframe+1,3)=length(find(OvertakeTip==j+0.1));
    end
    TakeOverTip=find(smooth(TailOut(:,3))>smooth(TailOut(:,2))); 
    TakeOverTail=find(smooth(TailOut(:,2))>smooth(TailOut(:,3)));
       
    if length(TakeOverTip)>length(TakeOverTail)
        if length(find(TailOut(:,2)>TailOut(:,3)))>max([0.05*length(TakeOverTip) 10])
            FilEnd(i).FinalImage=ImageTail;
            FilEnd(i).reject=1;
        else
            ImageUse=ImageTail;
            FilEnd(i).FinalImage=ImageTail;
        end
    else
        if length(find(TailOut(:,3)>TailOut(:,2)))>max([0.05*length(TakeOverTail) 10])
            FilEnd(i).FinalImage=ImageTip;
            FilEnd(i).reject=1;
        else
            ImageUse=ImageTip;
            FilEnd(i).FinalImage=ImageTip;
        end
        
    end
end

