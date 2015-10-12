%% convert mat file to netcdf
%----------------------------
clear all;close all;

matfile_dir  = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
matfile_name = 'ARISEairprocessed_with_insitu_withWVparams20150318_w_anacompare_w_consolidatedclouds20150318.mat';
matfile      = load([matfile_dir matfile_name]);

% create a combined struct to convert to netcdf

met.prof_20140904_1 = matfile.met20140904_.profnum1;
met.prof_20140904_2 = matfile.met20140904_.profnum2;
met.prof_20140904_3 = matfile.met20140904_.profnum3;
met.prof_20140904_4 = matfile.met20140904_.profnum4;
met.prof_20140904_5 = matfile.met20140904_.profnum5;
met.prof_20140904_6 = matfile.met20140904_.profnum6;

met.prof_20140906_1 = matfile.met20140906_.profnum1;
met.prof_20140906_2 = matfile.met20140906_.profnum2;

met.prof_20140907_1 = matfile.met20140907_.profnum1;
met.prof_20140907_2 = matfile.met20140907_.profnum2;
met.prof_20140907_3 = matfile.met20140907_.profnum3;
met.prof_20140907_4 = matfile.met20140907_.profnum4;
met.prof_20140907_5 = matfile.met20140907_.profnum5;
met.prof_20140907_6 = matfile.met20140907_.profnum6;

% create corresponding nc filename

met_nc_name = strcat(strtok(matfile_name,'.'),'.nc');

% convert moc to netcdf file

struct2nc(met,met_nc_name,'classic',0);

% open saved nc file

ncid    = netcdf.open(met_nc_name);
fncinfo = ncinfo(met_nc_name);

% theResult = mat2nc(theMatFile, theNetCDFFile, uniqueDims, noSqueeze)
% uniqueDims is 0 default; 1 means that each variable will get unique dims
% noSqueeze is 0 as default; 1 means that all singelton dimensions are
% remained intact
% this one has issues with mex files (too old)
% moc_nc = mat2nc(moc_combined_name, moc_nc_name);

