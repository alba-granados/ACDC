%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L1B/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L1B products:
% isardSAT and ESA for CryoSat-2 data runs the validation and saves the
% results into a .mat with a given structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_L1B_validation(filesBulk,i_files_valid,input_path_L1B_ISD,...
                                  input_path_L1B_ESA,input_path_L2_ESA,...
                                  path_comparison_results,surface_aligned,retrack_flag,varargin)

%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<8 || nargin>(8+2*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('input_path_L2_STL',{''},@(x)ischar(x));
p.addParamValue('figures_visible',0);
p.parse(varargin{:});
input_path_L2_STL=char(p.Results.input_path_L2_STL);
figures_visible=(p.Results.figures_visible);
clear p;
                              
filename_L1B_ISD=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
data_string=filename_L1B_ISD(17:17+30);
filename_L1B_ISD=strcat(input_path_L1B_ISD,filename_L1B_ISD);

inputL1BESAFiles   = dir(fullfile(input_path_L1B_ESA,['*' data_string(1:15) '*.DBL']));
inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' data_string(1:15) '*.DBL']));
if isempty(input_path_L2_STL)
    inputL2STLFiles   = 1;
else
    inputL2STLFiles   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
end

indexL1BESA   =   not([inputL1BESAFiles.isdir]);
filename_L1B_ESA=strcat(input_path_L1B_ESA,char(inputL1BESAFiles(indexL1BESA).name));
disp(char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name))
indexL2ESA   =   not([inputL2ESAFiles.isdir]);
filename_L2_ESA=strcat(input_path_L2_ESA,char(inputL2ESAFiles(indexL2ESA).name));
if ~isempty(input_path_L2_STL)
    indexL2STL   =   not([inputL2STLFiles.isdir]);
    filename_L2_STL=strcat(input_path_L2_STL,char(inputL2STLFiles(indexL2STL).name));
    [~]=L1B_validation_EM(filename_L1B_ISD,filename_L1B_ESA,filename_L2_ESA,path_comparison_results,surface_aligned,retrack_flag,...
                          'filename_L2_STL',filename_L2_STL,'figures_visible',figures_visible);
else
    [~]=L1B_validation_EM(filename_L1B_ISD,filename_L1B_ESA,filename_L2_ESA,path_comparison_results,surface_aligned,retrack_flag,...
                           'figures_visible',figures_visible);
end

end