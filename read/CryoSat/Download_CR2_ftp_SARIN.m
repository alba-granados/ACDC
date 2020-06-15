% BULK CR2 in your ass a FULL
% Download FBR and L1B products from a given list of L2 files or whatever.
% Inputs type of data needed and type of data provided
% download_CR2(downloaded,input)
% ftp= 'SIR_SAR_FR/2015/01/';

%%
% Important note: The FTP connections DOES NOT WORK with the VPN activated.

%%
%  close(mw);
ftpPath     = 'science-pds.cryosat.esa.int';
user        = 'cryosat081';
pswrd       = '4MsIuLNr';
downloading = 'SIR_SIN_L1';
inputDir    = 'F:\MoGLA\Central_Asia\L2_HDR\2015\';
outputDir   = 'F:\MoGLA\Central_Asia\inputs_2015\done\';
cd(outputDir);
mw = ftp(ftpPath,user,pswrd);


inputFiles      =   dir([inputDir '*.HDR']);
indexFiles      =   find(not([inputFiles.isdir]));
nFiles          =   length(indexFiles);


for iFile=1:nFiles
% for iFile=471:nFiles

%     progressbar(iFile/nFiles);
    input_filename = [inputFiles(indexFiles(iFile)).name];
    year    = input_filename(20:23);
    month   = input_filename(24:25);
    day     = input_filename(26:27);
    timeT1  = input_filename(29:34);
    cd(mw,['/' downloading '/' year '/' month '/']);
    try
        inputftpFiles   = dir(mw,['/' downloading '/' year '/' month '/*' year month day 'T' timeT1(1:end-2) '*.DBL']);
    catch
        try % 2nd try
            inputftpFiles   = dir(mw,['/' downloading '/' year '/' month '/*' year month day 'T' timeT1(1:end-2) '*.DBL']);
        catch
            disp(['Error dir: /' downloading '/' year '/' month '/*' year month day 'T' timeT1(1:end-2) '*.DBL'])
            continue;
        end
    end
    if(isempty(inputftpFiles))
        disp(['Not found: /' inputFiles(indexFiles(iFile)).name])
        continue;
    else
        indexftpFiles   =   find(not([inputftpFiles.isdir]));           
        nftpFiles       =   length(indexftpFiles);
        try
            mget(mw, inputftpFiles(indexftpFiles).name);    
        catch
            try % 2nd try
                mget(mw, inputftpFiles(indexftpFiles).name);
            catch
                disp(['Error downloading: /' inputFiles(indexFiles(iFile)).name])
                continue;
            end
        end
         
%          untar(inputftpFiles(indexftpFiles).name);
%          delete(inputftpFiles(indexftpFiles).name);
    end
    
    
end

