%% Motility of Filament Ends software
% user can select single or multiple files for analysis. parameter settings
% are set in MFE_Process.m

clearvars
close all

%% setting parameters

tic
%% end setting parameters

% Loading video file
[DirFile]=uipickfiles('filterspec','L:\*.avi');
for i=1:length(DirFile)
    dummy=strfind(DirFile{i},'\');
    filelist{i}=DirFile{i}(dummy(end)+1:end);
    dirout{i}=DirFile{i}(1:dummy(end));
end

if length(filelist)>1
    %iscell(filelist)
%     if ~isempty(dirout)
        parfor VideoNow=1:length(filelist)
            tic
            
            try
                MFE_Process(dirout{VideoNow},filelist{VideoNow},VideoNow)
            catch
                ['error in ' dir file]
            end
            toc
        end
%     end
else
    file=filelist;
    try
        MFE_Process(dirout{1},filelist{1},1)
    catch
        ['error in ' filelist{1}]
    end
    
end
toc