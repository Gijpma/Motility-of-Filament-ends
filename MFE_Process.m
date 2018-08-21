function MFE_Process(dirout,file,VideoNow)
%% setting parameters

param.MaxIter=0.4;          % maximum fraction of trace that can be a branch (default 0.2)
param.MinPresence=30;       % Remove filaments that are present for less than x frames
param.MaxGap=3;             % maximum number of frames an associated filament end is not present/detected
param.MaxDistAllow=8;       % maximum allowed filament end distance (default 8)
                            % decreasing this decreases the number of viable
                            % filaments, but increases the total number of traces
param.Square=14;            % size of square frame section around filament ends
param.BrightAdjust=.05;     %fraction adjustment for brightness cutoff
param.VelocityFactor=0.124*30;  % factor correcting pixel movement to velocity. scale
                            % is 0.124 mum per pixel, 30 frames per second
                            % video.
param.errormax=2;           % determines the closeness of fit of the segmentation algorithm. usually 1-2
param.CutOffVelocity=0.03;  % in mum/s   .03 for tonic, .15 testfiles
param.smooth=20;            % smoothing range, default 20 frames
param.MinSize=5;            % minimum size filament, in pixels. default 5
param.FrameAverage=2;       % number of frames to average. increase if noise dominates signal. default 2
param.DefaultThreshold=40;  % threshold value if no noise is present. default 40
param.MaxLength=9000;       % maximum duration for trace. low increases successful filament count
                            % but reduces ability to detect stopped vs. go,
                            % especially at low velocities.
%% end setting parameters

[frames, FilEnd, SizeFrame, num, Threshold,A]=VideoProcess(dirout, file, param);

FilEnd=TraceImage(FilEnd, SizeFrame);

for i=1:length(FilEnd)
    if FilEnd(i).reject==0
        [FilEnd(i).Frames, FilEnd(i).Dist,FilEnd(i).Trace, KillArray(i)]=TraceCreate(FilEnd(i).FinalImage,param.MaxIter, FilEnd(i).firstframe,FilEnd(i).lastframe, param);
        if KillArray(i)==1;
            FilEnd(i).reject=2;
        end
    end
end

save([dirout '\results\' file(1:end-4) '_ANL_Traced.mat'],'FilEnd','param','dirout','file','frames','VideoNow','num','Threshold','A')
%% end clean up trace

[TraceOut, FinalSpeeds, FinalPosx, FinalPosy, medianlength, FinalLength, Lengths]=TracePostProcess(FilEnd, param, frames);

[Velocity, fmot, fnumber]=SumVelocity(FinalSpeeds);

figure(VideoNow)
hold on
plot(smooth(Velocity,2),'b','displayname',file)
plot(smooth(fmot(5:end-5),3),'r')
saveas(gca,[dirout '\results\' file(1:end-4) '_result.png'])
saveas(gca,[dirout '\results\' file(1:end-4) '_result.fig'])
figure(VideoNow+100)
for i=1:length(FinalSpeeds(1,:))
    Moves(i)=length(find(FinalSpeeds(:,i)>0));
    Stops(i)=length(find(FinalSpeeds(:,i)<0));
end
plot(Moves,'g')
hold on
plot(Stops,'r')
saveas(gca,[dirout '\results\' file(1:end-4) '_result_filno.png'])
saveas(gca,[dirout '\results\' file(1:end-4) '_result_filno.fig'])
close(VideoNow)
close(VideoNow+100)
save([dirout '\results\' file(1:end-4) '_ANLpost.mat'],'TraceOut','FinalSpeeds','dirout','file','frames','FinalPosx','FinalPosy','Velocity','fmot')

VideoOut(FinalPosx, FinalPosy, FinalSpeeds, dirout, file);