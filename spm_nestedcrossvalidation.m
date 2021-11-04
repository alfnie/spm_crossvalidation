function spm_nestedcrossvalidation(SPM,xSPM,filename)
% SPM_NESTEDCROSSVALIDATION performs leave-one-out nested cross-validation of SPM second-level results
%
% When viewing second-level results in the standard spm_results window,
% type spm_nestedcrossvalidation; in the command line. Then select a .nii file
% where the cross-validation masks will be saved. The current second-level
% analyses and contrast will be re-evaluated leaving two subjects out at a
% time. For each of these cases supra-threshold voxels will be saved in a
% subject-specific mask file. 
%
% Description:
%  Within each (outer-layer) crossvalidation iteration the matrix data (of size number-of-subjects minus one by number-of-regions) will contain 
%  the cross-validated data of the selected second-level analysis (SPM.xY volumes) aggregated across each ROI. A linear model will be
%  used to fit this crossvalidated data to the design contrast values (X*C) for these same N-1 subjects, and this fit will be evaluated in the 
%  (outer-layer) left-out subject. The procedure will be repeated for each (outer-layer) iteration and the fit RMS and R2 values are reported
%
% Equations:
%   (1) Original SPM model:
%     Y = X*B + noise           where Y = functional data at each voxel; X = design matrix; B = estimated regressor coefficient images
%   (2) Hypothesis testing:
%     C*B = 0                   where C = between-subjects contrast
%   (3) Cross-validated data
%     y = Y*mask                where mask is a vector of 0/1 values indicating suprathreshold voxels where the above hypothesis was rejected (optionally a voxels x ROIs matrix, when dividing into multiple ROI parcels)
%   (4) Predictive model:
%     X*C' = a + y*b + noise    where a,b = estimated regressor coefficient values
%
% Procedure:
%   1) Using N-2 subjects (leaving subject i and j out), fit equations (1)&(2) above and compute mask_(i,j) (suprathreshold mask of voxels defined from an original SPM model that included all subjects except the i-th and j-th subjects)
%        note: for each j, mask_(i~=j,j) represents the masks that would result from running a standard spm_crossvalidation step on a restricted dataset that does not include subject #j
%                          mask_(j,j) represents the mask that would result from running the original SPM model on a restricted dataset that does not include subject #j
%   2) for each j,
%           2.1) for each i, compute equation (3) to estimate y_(i,j) = Y(i,:)*mask_(i,j)
%           2.2) use equation (4) restricted to all i~=j subjects to estimate a_j and b_j using linear regression (i.e. fit a N-1 subjects model of the form "X(i,:)*C' = a_j + y_(i,j)*b_j + noise" including all subjects i~=j)
%           2.3) compute error term E_j = a_j + y_(j,j)*b_j - X(j,:)*C on the left-out subject j
%   3) aggregate across all error terms to compute
%           rms = SUM_j=1:N{E_j^2} / N
%       
%   Step 1 above is performed by spm_nestedcrossvalidation.m
%   Steps 2 and 3 above are performed by spm_nestedcrossvalidation_extract (as equation (3) cross-validated data extraction can vary depending on choice of ROIs, a priori coordinates, etc.)
%
% see also SPM_NESTEDCROSSVALIDATION_EXTRACT
%

% alfnie@gmail.com 2021

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
V=repmat(V,[N^2,1]);for n=1:N^2,V(n).n=[n,1];end
V=spm_create_vol(V);
XYZmm=cell(1,N^2);
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
for N1=1:N
    for N2=1:N
        n=N1+N*(N2-1);
        cd(tempdir);
        if DOFAST % fast version
            if DOMULTIVARIATE
                X=[SPM.xX_multivariate.X,zeros(N,2)];
                X(N1,end-1)=1;
                X(N2,end)=1;
                if u<1,
                    [h,f,p,dof,stats]=conn_glm(X,Y,[SPM.xX_multivariate.C zeros(size(SPM.xX_multivariate.C,1),2)],SPM.xX_multivariate.M);
                    idx=find(p<u);
                else
                    [h,f]=conn_glm(X,Y,[SPM.xX_multivariate.C zeros(size(SPM.xX_multivariate.C,1),2)],SPM.xX_multivariate.M);
                    idx=find(f>u);
                end
            else
                X=[SPM.xX.X,zeros(N,2)];
                X(N1,end-1)=1;
                X(N2,end)=1;
                if u<1,
                    [h,f,p,dof,stats]=conn_glm(X,Y(:,:),[SPM.xCon(xSPM.Ic).c' 0 0],[],'collapse_none');
                    idx=find(p<u);
                else
                    [h,f]=conn_glm(X,Y(:,:),[SPM.xCon(xSPM.Ic).c' 0 0],[],'collapse_none');
                    idx=find(f>u);
                end
            end
            M=zeros(a(1).dim);
            M(idx)=1;
            cd(filepath); spm_write_vol(V(n),M);
            axes(hax);imagesc(mean(M,3));axis equal tight; title(hax,['Mask for subject #',num2str(N1),'&',num2str(N2),' (',num2str(nnz(M)),' voxels)']);axis equal;axis off;drawnow;
            fprintf('Subject %d&%d: %d supra-threshold voxels\n',N1,N2,nnz(M));
        else % slow version
            error('not implemented');
        end
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
fprintf('Subject-specific nested cross-validation mask file saved as %s\n',fullfile(filepath,[filename,'.nii']));

% extract cross-validated data
spm_nestedcrossvalidation_extract(fullfile(filepath,[filename,'.nii']),'?',fullfile(SPM.swd,'SPM.mat'));
