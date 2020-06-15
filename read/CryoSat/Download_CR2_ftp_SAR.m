% BULK CR2 in your ass a FULL
% Download FBR and L1B products from a given list of L2 files or whatever.
% Inputs type of data needed and type of data provided
% download_CR2(downloaded,input)
% ftp= 'SIR_SAR_FR/2015/01/';
% close(mw);
ftpPath     = 'science-pds.cryosat.esa.int';
user        = 'cryosat081';
pswrd       = '4MsIuLNr';
downloading = 'SIR_SAR_FR';
inputDir    = 'C:\Users\albert\Documents\WORK\SEOM\SHAPE\DATA\SAR\Vanern\L2_hdr\';
outputDir   = 'C:\Users\albert\Documents\WORK\SEOM\SHAPE\DATA\SAR\Vanern\FBR\';
cd(outputDir);
mw = ftp(ftpPath,user,pswrd);

FILES_COMPRESSED = 1;

inputFiles      =   dir([inputDir '*.HDR']);
indexFiles      =   find(not([inputFiles.isdir]));
nFiles          =   length(indexFiles);


for iFile=1:nFiles

    progressbar(iFile/nFiles);
    input_filename = [inputDir inputFiles(indexFiles(iFile)).name];
    input_filename = [inputFiles(indexFiles(iFile)).name];
    year    = input_filename(20:23);
    month   = input_filename(24:25);
    day     = input_filename(26:27);
    timeT1  = input_filename(29:34);
    cd(mw,['/' downloading '/' year '/' month '/']);
    
    if FILES_COMPRESSED == 1
        inputftpFiles   = dir(mw,['/' downloading '/' year '/' month '/*' year month day 'T' timeT1 '*']);
    else
        inputftpFiles   = dir(mw,['/' downloading '/' year '/' month '/*' year month day 'T' timeT1 '*.DBL']);
    end
    
    %%
    if(~isempty(inputftpFiles))
        indexftpFiles   =   find(not([inputftpFiles.isdir]));           
        nftpFiles       =   length(indexftpFiles);
         mget(mw, inputftpFiles(indexftpFiles).name,outputDir);
         if FILES_COMPRESSED == 1
             untar(inputftpFiles(indexftpFiles).name);
             delete(inputftpFiles(indexftpFiles).name);
         end
    end
    
    
end

