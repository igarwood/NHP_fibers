function [session_info] = setup(session)


function_file = mfilename('fullpath');
[path,~,~] = fileparts(function_file);
dir_end = strfind(path,'/');
directory = path(1:dir_end(end));
folder = [directory,'experimental_data/'];
device_folder = [directory,'device_data/'];
save_folder = [directory,'processed_data/'];

fs_lfp = 1000;
fs_spike = 30000;

session_info.folder = folder;
session_info.device_folder = device_folder;
session_info.save_folder = save_folder;
session_info.fs_lfp = fs_lfp;
session_info.fs_spike = fs_spike;
if nargin == 1
    session_info.fileroot = session;
    session_info.pharmacology = extract_pharmacology([folder,session]);
end
end
