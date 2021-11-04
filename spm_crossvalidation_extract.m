function [data,voxels,contrasts,labels]=spm_crossvalidation_extract(maskname,roiname,spmname)
% SPM_CROSSVALIDATION_EXTRACT extracts BOLD data from subject-specific
% cross-validation mask file. 
%
% After running SPM_CROSSVALIDATION use
%    data=SPM_CROSSVALIDATION_EXTRACT;
% to extract cross-validated data using the resulting subject-specific masks.
%
% Steps:
%    1) select a cross-validation mask file (*.nii file created using SPM_CROSSVALIDATION)
%    2) select ROI file(s) (this file(s) should contain a priori ROI definitions,
%    the data from each subject-specific mask will be aggregated across
%    each ROI; you can leave this empty if you prefer to aggregated across
%    ALL of the voxels in the subject-specific masks)
%    3) select second-level SPM.mat file (typically the same second-level
%    SPM.mat file used in the SMP_CROSSVALIDATION step)
%
% The matrix data (of size number-of-subjects by number-of-regions) will contain 
% the cross-validated data of the selected second-level analysis (SPM.xY volumes) 
% aggregated across each ROI. 
%
% see also SPM_CROSSVALIDATION
%

% alfnie@gmail.com 06/10

if nargin<1
    maskname=spm_select(1,'\.nii$|\.img$','Select cross-validation mask file');
end
[maskpath,maskname,maskext]=fileparts(maskname);
if nargin<2||(ischar(roiname)&&strcmp(roiname,'?')),
    toptions={'All suprathreshold voxels','Largest N clusters','All voxels within user-defined ROI file','Closest cluster to XYZ coordinates'};
    [idx,ok]=listdlg('PromptString','Extract from:','ListString',toptions,'SelectionMode','single','ListSize',[300 100]);
    if ~ok, return; end
    switch(idx)
        case 1, roiname=[];
        case 2, roiname=inputdlg({'Number of clusters'},'',1,{'1'}); roiname=str2num(roiname{1}); if numel(roiname)~=1, error('incorrect value'); end
        case 3, roiname=spm_select(inf,'\.tal$|\.img$|\.nii$','Select ROI file(s)');
        case 4, roiname=inputdlg({'X Y Z coordinates (mm)'},'',1,{'0 0 0'}); roiname=str2num(roiname{1}); if numel(roiname)~=3, error('incorrect value'); end
    end
end
if nargin<3
    if ~isempty(dir(fullfile(maskpath,[maskname,'.cv']))),
        load(fullfile(maskpath,[maskname,'.cv']),'SPM','-mat');
        spmname={fullfile(SPM.swd,'SPM.mat')};
    else spmname='';end
    spmname=spm_select(1,'SPM.*\.mat$','Select second-level analysis file',spmname);
end

%load(fullfile(maskpath,[maskname,'.cv']),'SPM','-mat');
load(spmname,'SPM');
N=size(SPM.xX.X,1);
ni=nifti([maskname maskext]);
Nm=[ni.dat.dim 1 1 1 1];
Nm=Nm(4);

if ~isempty(roiname)
    Vextract=struct('fname',fullfile(maskpath,[maskname,'.extracted',maskext]),'descrip','spm_crossvalidation_extract file','mat',ni.mat,'dim',ni.dat.dim(1:3),'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('uint8'),spm_platform('bigend')]);
    try, spm_unlink(Vextract.fname); end
    Vextract=repmat(Vextract,[N,1]);for nb=1:N,Vextract(nb).n=[nb,1];end
    Vextract=spm_create_vol(Vextract);
end

data=[];voxels=[];
if ~isfield(SPM.xY,'P')
    SPM.xY.P=cellfun(@(a,b)[a,',',b],{SPM.xY.VY.fname},cellfun(@(x)num2str(x(1)),{SPM.xY.VY.n},'uni',0),'uni',0);
end
if isempty(roiname) % extract average across ALL suprathrehsold voxels
    for n=1:N,
        m=1+mod(n-1,Nm);
        disp(['Extracting data from subject ',num2str(m)]);
        [data(n,:),labels]=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'level','rois');
        voxels(n,:)=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'summary_measure','size','level','rois');
    end
elseif ischar(roiname) % extract average across all suprathreshold voxels within each ROI defined in the file 'roiname'
    for n=1:N,
        m=1+mod(n-1,Nm);
        disp(['Extracting data from subject ',num2str(m)]);
        [data(n,:),labels,tparams]=rex(SPM.xY.P{n},roiname,'level','clusters','conjunction_mask',fullfile(maskpath,[maskname,maskext,',',num2str(m)]));
        voxels(n,:)=rex(SPM.xY.P{n},roiname,'summary_measure','size','level','clusters','conjunction_mask',fullfile(maskpath,[maskname,maskext,',',num2str(m)]));
        xyz=cat(1,tparams.ROIinfo.voxels{1}{:}); 
        temp=zeros(ni.dat.dim(1:3)); temp(1+(round([xyz ones(size(xyz,1),1)]*pinv(ni.mat)')-1)*[1;ni.dat.dim(1);ni.dat.dim(1)*ni.dat.dim(2);0])=1; Vextract(n)=spm_write_vol(Vextract(n),temp); 
    end
elseif isnumeric(roiname)&&numel(roiname)==1 % extract average across all suprathreshold voxels within each of the N largest clusters
    for n=1:N,
        m=1+mod(n-1,Nm);
        disp(['Extracting data from subject ',num2str(m)]);
        [tdata{n},labels,tparams]=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'level','clusters');
        tvoxels{n}=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'summary_measure','size','level','clusters');
        xyz=cat(1,tparams.ROIinfo.voxels{1}{1:min(numel(tparams.ROIinfo.voxels{1},roiname))}); 
        temp=zeros(ni.dat.dim(1:3)); temp(1+(round([xyz ones(size(xyz,1),1)]*pinv(ni.mat)')-1)*[1;ni.dat.dim(1);ni.dat.dim(1)*ni.dat.dim(2);0])=1; Vextract(n)=spm_write_vol(Vextract(n),temp); 
    end
    data=cellfun(@(x)[x(1:min(numel(x),roiname)) nan(1,max(0,roiname-numel(x)))],tdata,'uni',0); % selects N ('roiname') largest clusters
    data=cell2mat(data(:));
    voxels=cellfun(@(x)[x(1:min(numel(x),roiname)) nan(1,max(0,roiname-numel(x)))],tvoxels,'uni',0);
    voxels=cell2mat(voxels(:));
elseif isnumeric(roiname)&&numel(roiname)==3 % extract average across cluster closest to XYZ coordinates
    for n=1:N,
        m=1+mod(n-1,Nm);
        disp(['Extracting data from subject ',num2str(m)]);
        [tdata{n},labels,tparams]=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'level','clusters');
        tvoxels{n}=rex(SPM.xY.P{n},fullfile(maskpath,[maskname,maskext,',',num2str(m)]),'summary_measure','size','level','clusters');
        tdist=[]; for n1=1:numel(tparams.ROIinfo.voxels{1}), tdist(n1)=min(sum(abs(tparams.ROIinfo.voxels{1}{n1}-repmat(roiname,size(tparams.ROIinfo.voxels{1}{n1},1),1)).^2,2)); end
        [nill,idx]=sort(tdist);
        tdata{n}=tdata{n}(idx);
        tvoxels{n}=tvoxels{n}(idx);
        xyz=tparams.ROIinfo.voxels{1}{idx}; 
        temp=zeros(ni.dat.dim(1:3)); temp(1+(round([xyz ones(size(xyz,1),1)]*pinv(ni.mat)')-1)*[1;ni.dat.dim(1);ni.dat.dim(1)*ni.dat.dim(2);0])=1; Vextract(n)=spm_write_vol(Vextract(n),temp); 
    end
    data=cellfun(@(x)[x(1:min(numel(x),1)) nan(1,max(0,1-numel(x)))],tdata,'uni',0); % selects closest cluster
    data=cell2mat(data(:));
    voxels=cellfun(@(x)[x(1:min(numel(x),1)) nan(1,max(0,1-numel(x)))],tvoxels,'uni',0);
    voxels=cell2mat(voxels(:));
else
    error('incorrect field roiname (third input argument should be empty, a number, or a filename)');
end

KWY   = spm_filter(SPM.xX.K,SPM.xX.W*data);
beta  = SPM.xX.pKX*KWY;
clear contrasts;
for ncon=1:numel(SPM.xCon),
    y=SPM.xX.X*SPM.xCon(ncon).c; % estimate contrast sum for each subject
    x=[ones(N,1) data];          % linearly from the data
    test_y=nan(N,size(y,2));
    null_y=nan(N,size(y,2));
    for n=1:N,
        train_y=y([1:n-1,n+1:end],:);
        train_x=x([1:n-1,n+1:end],:);
        train_beta=train_x\train_y;
        null_beta=train_x(:,1)\train_y;
        test_y(n,:)=x(n,:)*train_beta;
        null_y(n,:)=x(n,1)*null_beta;
    end
    test_rms=sqrt(mean(mean(abs(y-test_y).^2)));
    null_rms=sqrt(mean(mean(abs(y-null_y).^2)));
    r2=max(0,1-test_rms^2/null_rms^2);
    contrasts(ncon)=struct('design',SPM.xX.X,'name',SPM.xCon(ncon).name,'STAT',SPM.xCon(ncon).STAT,'contrast',SPM.xCon(ncon).c,'values',SPM.xCon(ncon).c'*beta,'CV_Y',y,'CV_Yfit',test_y,'CV_RMS',test_rms,'CV_R2',r2);
    if size(data,2)==1, 
        disp(['Contrast ''',SPM.xCon(ncon).name,''' cross-validated value(s) : ',num2str(contrasts(ncon).values(:)')]); 
        disp(['Contrast ''',SPM.xCon(ncon).name,''' cross-validated linear fit R-square : ',num2str(contrasts(ncon).fit_R2(:)')]);
    end
end

if ~nargout
    assignin('base','CV_DATA',data);
    assignin('base','CV_CONTRASTS',contrasts);
    assignin('base','CV_VOXELS',voxels);
    assignin('base','CV_LABELS',labels);
    disp('Output stored in workspace''s variable CV_DATA, CV_CONTRASTS, CV_VOXELS, and CV_LABELS');
end
