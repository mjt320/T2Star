function pipeline_R2s_convert(opts);

mkdir(opts.niftiDir); delete([opts.niftiDir '/*.*']);

%% convert dicoms to 3D niftis
dicomPaths=getMultipleFilePaths([opts.dicomDir '/*.dcm']);
if isempty(dicomPaths); dicomPaths=getMultipleFilePaths([opts.dicomDir '/*.IMA']); end;
system(['dcm2niix -f magn_image_%e -t n -v n -b y -z n -o ' opts.niftiDir ' ' opts.dicomDir]);

%% make a 4D nifti containing only magnitude data
fileList={};
acqPars.NEchoes=size(dir([opts.niftiDir '/magn_image_*.nii']),1);
for iEcho=1:acqPars.NEchoes
    fileList=[fileList; opts.niftiDir '/magn_image_' num2str(iEcho) '.nii'];
end
spm_file_merge(fileList,[opts.niftiDir '/mag4D.nii'],0);

copyfile(fileList{1},[opts.niftiDir '/firstEcho.nii']); % copy 3D magnitude image of first echo
delete([opts.niftiDir '/magn_image_*.nii']); %delete 3D files

%% loop through dicoms and get echo times (dicominfo is surprisingly slow)
allTE=[];
while 1
    iDicom=floor(rand*size(dicomPaths,1))+1;
    temp=dicominfo(dicomPaths{iDicom}); allTE=[allTE; temp.EchoTime];
    if size(unique(allTE),1)==acqPars.NEchoes; break; end
end
acqPars.TE=0.001*sort(unique(allTE)); %get unique echo times (s) and store in ascending order
if size(acqPars.TE,1)~=acqPars.NEchoes; error('Number of echo times is not equal to expected number of echoes.'); end

save([opts.niftiDir '/acqPars'],'acqPars');

end