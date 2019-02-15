function VideoOut(FinalPosx, FinalPosy, FinalSpeeds, dirout, file)
% Generates an output video indicating moving (green) vs stopped (red)
% filaments and the velocity of the moving filaments as indicated by the
% hight of the green bar (scaled to the mean velocity in the video)
tic
VideoObject = VideoReader([dirout file]);
v = VideoWriter([dirout '\results2\' file(1:end-4) '_TipTail.avi'],'Motion JPEG AVI');
open(v)
vidFrames = readFrame(VideoObject);
factor=255/mean(FinalSpeeds(find(FinalSpeeds>0)));
FinalPosx(find(~(FinalPosx>-20)))=10;
FinalPosy(find(~(FinalPosy>-20)))=10;
for i=3:find(abs(sum(FinalSpeeds))>0,1,'last')-1        
    
    if mod(i,4)==0
        
        VideoObject.CurrentTime=i/30;
        vidFrames = readFrame(VideoObject);
        A=imresize(vidFrames,1);
        
        for j=1:find(abs(sum(FinalSpeeds'))>0,1,'last')        
            y=full(round(FinalPosy(j,i)));
            x=full(round(FinalPosx(j,i)));
            if FinalSpeeds(j,i)<0 & x>2 & y>2
                A(x-1:x+1,y-1:y+1,1)=255;
                A(x,y,2)=255;
                A(x,y,3)=255;
                
            elseif FinalSpeeds(j,i)>0 & FinalSpeeds(j,i)<50
                A(x-1:x+1,y-1:y+1,2)=full(round(FinalSpeeds(j,i)*factor));
                xout=x-1*full(round(FinalSpeeds(j,i)*factor/10));
                if xout<1
                    xout=1;
                end
                
                A(xout:x+1,y-1:y+1,2)=255;
                A(x,y,1)=255;
                A(x,y,3)=255;
%                 if j==16
%                     A(x-1:x+1,y-1:y+1,1)=0;
%                 end
            end
        end
       
%        A=imresize(A,[480 720]);
      
        writeVideo(v,A);
    end
    
end
close(v)
toc