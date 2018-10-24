function pipeline_R2s_create_map_new(opts)

load([opts.niftiDir '/acqPars']);

mkdir(opts.mapDir); delete([opts.mapDir '/*.*']);

if ~isfield(opts,'fit'); opts.fit=ones(1,acqPars.NEchoes); end %if opts.fit not specified, fit all data

%% load 4D magnitude data
[magnitude,xyz]=spm_read_vols(spm_vol([opts.niftiDir '/mag4D.nii']));

%% initialise output arrays
volTemplate=spm_vol([opts.niftiDir '/firstEcho.nii']); %use this header as template for 3D output files
R2s=nan(volTemplate.dim); S0=nan(volTemplate.dim); RSq=nan(volTemplate.dim); model=nan([volTemplate.dim sum(opts.fit,2)]);

%% do the fitting

for i3=1:size(magnitude,3); for i1=1:size(magnitude,1); for i2=1:size(magnitude,2) %loop through voxels
            
            [S0(i1,i2,i3),R2s(i1,i2,i3),T2s(i1,i2,i3),RSq(i1,i2,i3),model(i1,i2,i3,:)] = ...
                fit_R2s(acqPars.TE,magnitude(i1,i2,i3,:),opts.fit,opts.threshold1,opts.threshold2);
            
        end;
    end;
    disp([num2str(i3) '/' num2str(size(magnitude,3))]);
end;

%% write output images
SPMWrite4D(volTemplate,model,opts.mapDir,'model',16);
SPMWrite4D(volTemplate,magnitude(:,:,:,opts.fit==1),opts.mapDir,'signal',16);

paramNames={'S0' 'R2s' 'T2s' 'RSq'};
outputs={S0 R2s T2s RSq};

for iOutput=1:size(outputs,2)
    SPMWrite4D(volTemplate,outputs{iOutput},opts.mapDir,paramNames{iOutput},16);
end

end