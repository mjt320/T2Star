function pipeline_R2s_convert(opts);

mkdir(opts.niftiDir); delete([opts.niftiDir '/*.*']);

%% convert dicoms to 3D niftis
dicomPaths=getMultipleFilePaths([opts.dicomDir '/*.dcm']);
if isempty(dicomPaths); dicomPaths=getMultipleFilePaths([opts.dicomDir '/*.IMA']); end;

MLB{1}.spm.util.dicom.data = dicomPaths;
MLB{1}.spm.util.dicom.root = 'flat';
MLB{1}.spm.util.dicom.outdir = {opts.niftiDir};
MLB{1}.spm.util.dicom.convopts.format = 'nii';
MLB{1}.spm.util.dicom.convopts.icedims = 0;
spm_jobman('run', MLB);

%% make a 4D nifti containing only magnitude data
fileList={};
for iEcho=1:opts.NEchoes
    temp=sort(getMultipleFilePaths([opts.niftiDir '/*' num2str(iEcho,'%02d') '.nii'])); %get all files associated with this echo
    echoImageFile=temp{1}; %assume first file is magnitude 
    fileList=[fileList; echoImageFile];
end
spm_file_merge(fileList,[opts.niftiDir '/mag4D.nii'],0);

copyfile(fileList{1},[opts.niftiDir '/firstEcho.nii']); % copy 3D magnitude image of first echo
for iEcho=1:opts.NEchoes; delete([opts.niftiDir '/*-' num2str(iEcho,'%02d') '.nii']); end %delete 3D files

%% loop through dicoms and get echo times (dicominfo is surprisingly slow)
allTE=[];
while 1
    iDicom=floor(rand*size(dicomPaths,1))+1;
    temp=dicominfo(dicomPaths{iDicom}); allTE=[allTE; temp.EchoTime];
    if size(unique(allTE),1)==opts.NEchoes; break; end
end
acqPars.TE=0.001*sort(unique(allTE)); %get unique echo times (s) and store in ascending order
if size(acqPars.TE,1)~=opts.NEchoes; error('Number of echo times is not equal to expected number of echoes.'); end

save([opts.niftiDir '/acqPars'],'acqPars');

end