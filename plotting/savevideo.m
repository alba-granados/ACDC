% savevideo(workingDir)

%workingDir = '/net/n8800pro/share/roger/Jason-CS/Plots/GPP_Ph3/Scenario3/Videos/4_Ambiguity_Analysis';
% workingDir = './results/stacks2/';
workingDir = ['./'];

% mkdir(workingDir);
% mkdir(workingDir,'images');

outputVideo = VideoWriter(fullfile(workingDir,'TRP_rotation_HD.avi'));
outputVideo.FrameRate = 25;
open(outputVideo);

imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';

for i_imag = 1:length(imageNames)
    img = imread(fullfile(workingDir,imageNames{i_imag}));

    writeVideo(outputVideo,img);
end
close(outputVideo);