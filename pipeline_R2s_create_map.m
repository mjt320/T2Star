function pipeline_R2s_create_map(opts)

load([opts.niftiDir '/acqPars']);

mkdir(opts.mapDir); delete([opts.mapDir '/*.*']);

if ~isfield(opts,'fit'); opts.fit=ones(1,acqPars.NEchoes); end %if opts.fit not specified, fit all data

%% load 4D magnitude data
[magnitude,xyz]=spm_read_vols(spm_vol([opts.niftiDir '/mag4D.nii']));

%% initialise output arrays
volTemplate=spm_vol([opts.niftiDir '/firstEcho.nii']); %use this header as template for 3D output files
R2s=nan(volTemplate.dim); S0=nan(volTemplate.dim); RSq=nan(volTemplate.dim); model=nan([volTemplate.dim sum(opts.fit,2)]);

%% do the fitting
s=@(coef,t) coef(1)*exp(-coef(2)*t); % this is the signal model function

for i3=1:size(magnitude,3); for i1=1:size(magnitude,1); for i2=1:size(magnitude,2) %loop through voxels
            
            if magnitude(i1,i2,i3,1)<opts.threshold1; continue; end %skip low intensity voxels
            
            tic;
            t=acqPars.TE(opts.fit==1); %define independent variable (time)
            y=squeeze(magnitude(i1,i2,i3,opts.fit==1)); %define dependent variable (signal)
            
            %% use linear regression to determine initial guess            
            if size(y,1)>2; NLinReg=3; else NLinReg=2; end
            temp=regress(log(y(1:NLinReg)),[ones(NLinReg,1) t(1:NLinReg)]);
            x0=[exp(temp(1)) -temp(2)];
            if x0(1)<0 || imag(x0(1))~=0 || isnan(x0(1)) || isinf(x0(1))...
                    || x0(2)<0 || imag(x0(2))~=0 || isnan(x0(2)) || isinf(x0(2)); x0=[1000 10]; end
            
            
            %% use non-linear least squares to fit
            typicalX=[max(y) 50];
            
            % use only data points above opts.threshold2 to avoid fitting noise
            % if there are fewer than 2 data points selected then don't fit
            idx_reduced=y>opts.threshold2;
            if sum(idx_reduced)<2; continue; end
            
            [x,resnorm,residual,exitflag,output]=lsqcurvefit(s,x0,t(idx_reduced),y(idx_reduced)...
                ,[],[],optimset('Display','off','TypicalX',typicalX,'Algorithm','levenberg-marquardt'));
            
            S0(i1,i2,i3)=x(1); R2s(i1,i2,i3)=x(2);
            model(i1,i2,i3,idx_reduced)=s(x,t(idx_reduced));
            RSq(i1,i2,i3)=1 - sum((y(idx_reduced).'-squeeze(model(i1,i2,i3,idx_reduced)).').^2) / (sum((y(idx_reduced).'-mean(y(idx_reduced))).^2));
            
            timeElapsed=toc;
            
            if rand<0.001 %randomly plot data to check it's working
                figure(1),plot(t,y,'ko',t(idx_reduced),squeeze(model(i1,i2,i3,idx_reduced)),'b-')
                title({['Initial coefficients: ' num2str(x0)] ['Fitted coefficients: ' num2str(x)] ['Time elapsed: ' num2str(timeElapsed)] ['Exit flag: ' num2str(exitflag)] ['Evaluations: ' num2str(output.funcCount)]});
                pause(0.1);
            end
        end;
    end;
    disp([num2str(i3) '/' num2str(size(magnitude,3))]);
end;

T2s=1./R2s;

%% write output images
SPMWrite4D(volTemplate,model,opts.mapDir,'model',16);
SPMWrite4D(volTemplate,magnitude(:,:,:,opts.fit==1),opts.mapDir,'signal',16);

paramNames={'S0' 'R2s' 'T2s' 'RSq'};
outputs={S0 R2s T2s RSq};

for iOutput=1:size(outputs,2)
    SPMWrite4D(volTemplate,outputs{iOutput},opts.mapDir,paramNames{iOutput},16);
end

end