function spm_crossvalidation(SPM,xSPM,filename)
% SPM_CROSSVALIDATION performs leave-one-out cross-validation of SPM second-level results
%
% When viewing second-level results in the standard spm_results window,
% type spm_crossvalidation; in the command line. Then select a .nii file
% where the cross-validation masks will be saved. The current second-level
% analyses and contrast will be re-evaluated leaving one subject out at a
% time. For each of these cases supra-threshold voxels will be saved in a
% subject-specific mask file. 
%
% see also SPM_CROSSVALIDATION_EXTRACT
%

% alfnie@gmail.com 06/10

DOFAST=true;
ok=1;
try
    if nargin<1||isempty(SPM), SPM=[]; SPM=evalin('base','SPM'); end
    if nargin<2||isempty(xSPM), xSPM=[]; xSPM=evalin('base','xSPM'); end
catch
    ok=0;
end
if ~ok||isempty(SPM)||isempty(xSPM), [SPM,xSPM]=spm_getSPM; end
if nargin<3||isempty(filename), 
    [filename,filepath]=uiputfile('*.nii','Save cross-validation mask file as','cvmask.nii'); 
    filename=fullfile(filepath,filename);
end

[filepath,filename,fileext]=fileparts(filename);
if isempty(filepath), filepath=pwd; end

try
[t1,t2,t3,t4]=strread(xSPM.thresDesc,'%c%c%f %s');
catch
    disp('tal');
    disp(ok);
    whos
    disp(evalin('base','SPM'));
    disp(SPM);
    disp(xSPM);
    disp(xSPM.thresDesc);
t=textscan(xSPM.thresDesc,'%c%c%f %s');
[t1,t2,t3,t4]=deal(t{:});
dbstack -completenames
end

if isempty(t4)&&all(t2=='='),u=t3; thresDesc='none';
else 
    switch(t4{1})
        case '(unc.)', u=t3; thresDesc='none';
        case '(FWE)',  u=t3; thresDesc='FWE';
        case '(FDR)',  u=t3; thresDesc='FDR';
        otherwise,     error(['unrecognized xSPM threshold format ',xSPM.thresDesc]);
    end
end
if isfield(SPM,'xX_multivariate')&&isfield(SPM.xX_multivariate,'M')&&~isequal(SPM.xX_multivariate.M,1), DOMULTIVARIATE=true;
else DOMULTIVARIATE=false;
end

cwd=pwd;
tempdir=fullfile(filepath,[filename,'_',char(64+ceil(6*rand(1,10)))]);
[ok,msg]=mkdir(tempdir); if ~ok||~isempty(msg), error(['unable to create directory ',tempdir]); end

% evaluate leave-one-out cross-validated models
if DOMULTIVARIATE, N=size(SPM.xX_multivariate.X,1);
else N=size(SPM.xX.X,1);
end
cd(filepath);
V=struct('fname',[filename,'.nii'],...
    'mat',SPM.xCon(1).Vcon.mat,...
    'dim',SPM.xCon(1).Vcon.dim,...
    'n',[1,1],...
    'pinfo',[1;0;0],...
    'dt',[spm_type('uint8') spm_platform('bigend')],...
    'descrip',sprintf('spm_crossvalidation (subject-specific crossvalidation mask, analysis %s',SPM.swd));
V=repmat(V,[N,1]);for n=1:N,V(n).n=[n,1];end
V=spm_create_vol(V);
XYZmm=cell(1,N);
I=eye(N);
Finter = spm_figure('FindWin','Interactive'); if ~isempty(Finter),set(Finter,'tag','temporal');end
figure('name',mfilename,'numbertitle','off','color','w','colormap',gray);hax=gca;axis off;
if ~(DOFAST&&strcmp(thresDesc,'none')&&numel(xSPM.Ic)==1&&(size(SPM.xCon(xSPM.Ic).c,2)==1||DOMULTIVARIATE)),DOFAST=false; end
if DOFAST
    if isfield(SPM.xY,'P'), a=spm_vol(char(SPM.xY.P));
    else a=SPM.xY.VY;
    end
    Y=permute(spm_read_vols(a),[4 1 2 3]);
    if DOMULTIVARIATE
        Y=reshape(Y,N,[],size(Y,2),size(Y,3),size(Y,4));
    end
end
for n=1:N
    cd(tempdir);
    if DOFAST % fast version
        if DOMULTIVARIATE
            X=[SPM.xX_multivariate.X,zeros(N,1)];
            X(n,end)=1;
            if u<1,
                [h,f,p,dof,stats]=conn_glm(X,Y,[SPM.xX_multivariate.C zeros(size(SPM.xX_multivariate.C,1),1)],SPM.xX_multivariate.M);
                idx=find(p<u);
            else
                [h,f]=conn_glm(X,Y,[SPM.xX_multivariate.C zeros(size(SPM.xX_multivariate.C,1),1)],SPM.xX_multivariate.M);
                idx=find(f>u);
            end
        else
            X=[SPM.xX.X,zeros(N,1)];
            X(n,end)=1;
            if u<1,
                [h,f,p,dof,stats]=conn_glm(X,Y(:,:),[SPM.xCon(xSPM.Ic).c' 0],[],'collapse_none');
                idx=find(p<u);
            else
                [h,f]=conn_glm(X,Y(:,:),[SPM.xCon(xSPM.Ic).c' 0],[],'collapse_none');
                idx=find(f>u);
            end
        end
        M=zeros(a(1).dim);
        M(idx)=1;
        cd(filepath); spm_write_vol(V(n),M);
        axes(hax);imagesc(mean(M,3));axis equal tight; title(hax,['Mask for subject #',num2str(n),' (',num2str(nnz(M)),' voxels)']);axis equal;axis off;drawnow;
        fprintf('Subject %d: %d supra-threshold voxels\n',n,nnz(M)); 
    else % slow version
        % model estimation
        tSPM=SPM; tSPM.swd=tempdir;
        spm_unlink(fullfile(tempdir,'mask.img'));
        spm_unlink(fullfile(tempdir,'mask.nii'));
        if ~isfield(tSPM.xVi,'Vi') || (isfield(tSPM.xVi,'form')&&strcmp(tSPM.xVi.form,'i.i.d.')),tSPM=rmfield(tSPM,'xVi');
        else tSPM.xVi=rmfield(tSPM.xVi,'V'); tSPM.xX=rmfield(tSPM.xX,'W'); end
        tSPM.xX.X=cat(2,tSPM.xX.X,I(:,n));
        tSPM.xX.name{end+1}=sprintf('spm_crossvalidation subject #%d',n);
        tSPM=spm_spm(tSPM);
        % contrast estimation
        clear xCon; for n1=1:numel(xSPM.Ic), xCon(n1) = spm_FcUtil('Set',SPM.xCon(xSPM.Ic(n1)).name,SPM.xCon(xSPM.Ic(n1)).STAT,'c',cat(1,SPM.xCon(xSPM.Ic(n1)).c,zeros(1,size(SPM.xCon(xSPM.Ic(n1)).c,2))),tSPM.xX.xKXs); end
        for n1=1:numel(xSPM.Im), xCon(numel(xSPM.Ic)+n1) = spm_FcUtil('Set',SPM.xCon(xSPM.Im(n1)).name,SPM.xCon(xSPM.Im(n1)).STAT,'c',cat(1,SPM.xCon(xSPM.Im(n1)).c,zeros(1,size(SPM.xCon(xSPM.Im(n1)).c,2))),tSPM.xX.xKXs); end
        tSPM.xCon=xCon;
        tSPM=spm_contrasts(tSPM,1:numel(xSPM.Ic)+numel(xSPM.Im));
        % threshold contrast
        txSPM=xSPM; txSPM.swd=tempdir;
        txSPM.u=u; txSPM.thresDesc=thresDesc;
        txSPM.Ic=1:numel(xSPM.Ic);
        txSPM.Im=numel(xSPM.Ic)+(1:numel(xSPM.Im));
        [tSPM,txSPM]=spm_getSPM(txSPM);
        % create mask file
        XYZmm{n}=txSPM.XYZmm; if size(XYZmm{n},1)<4, XYZmm{n}(4,:)=1; end
        idx=sub2ind(txSPM.DIM(:)',txSPM.XYZ(1,:),txSPM.XYZ(2,:),txSPM.XYZ(3,:));
        M=zeros(txSPM.DIM(:)');
        M(idx)=1;
        cd(filepath); spm_write_vol(V(n),M);
        axes(hax);spm_mip(txSPM.Z,txSPM.XYZmm,txSPM.M);title(hax,['Mask for subject #',num2str(n),' (',num2str(size(XYZmm{n},2)),' voxels)']);axis equal;axis off;drawnow;
        fprintf('Subject %d: %d supra-threshold voxels\n',n,size(XYZmm{n},2));
    end
end
if ~isempty(Finter),set(Finter,'tag','Interactive');end
% remove files
cd(tempdir);
files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
         '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
         '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$','SPM.mat'};
for i=1:length(files)
    j = spm_select('List',tempdir,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end
% save results file
cd(cwd);
save(fullfile(filepath,[filename,'.cv']),'XYZmm','SPM','xSPM');
[nill,ok]=system(['rmdir ',tempdir]);
fprintf('Subject-specific cross-validation mask file saved as %s\n',fullfile(filepath,[filename,'.nii']));

% extract cross-validated data
spm_crossvalidation_extract(fullfile(filepath,[filename,'.nii']),'?',fullfile(SPM.swd,'SPM.mat'));
