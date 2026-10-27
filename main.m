%% New ST variogram and GIVD, GIVE for USA
% Each IGP variogram fitting update
%%read_filedate ('I:\主机材料_polyu PhD\Anisotropy\Anisotropy\USA and Canada_Analysis\sp3s','I:\主机材料_polyu PhD\Anisotropy\Anisotropy\USA and Canada_Analysis','sp3');

%% Step 1: Prepocessing: Detrending
%function final_model=st_variofit(SPF,STPF,n,) % n+1 epochs used in PF and variogram
disp('must consider the time range for data processing');% generally, GIVE and UIVE updata on 5 minutes (10 epoch)
%% Step 1: Preprocessing: IPPVD:m,r for iono model in spatial PF and spat-temp PF
% base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\Hongkong\download';
base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\USA Quiet2';
% 获取所有文件夹列表
folders = dir(base_path);
folders = folders([folders.isdir]); % 只保留目录
folders = folders(~ismember({folders.name}, {'.', '..'})); % 排除 '.' 和 '..'
parpool('local', 13);

K0 = 7; K1 = 30;
N  = K1 - K0 + 1;
ok     = false(N,1);
wsLog  = strings(N,1);
errLog = cell(N,1);      % 用 cell 存错误结构体最稳
parfor idx = 1:N
    k = K0 + idx - 1;
    workspace = fullfile(base_path, folders(k).name);
    wsLog(idx) = string(workspace);

    try
        outMat = fullfile(workspace, 'e_s_vD_ipp', 'vD_all.mat');
        if exist(outMat, 'file') == 2
            A = load(outMat);
            IPPVD = A.IPPVD;
        else
            IPPVD = extractIPPVD(workspace);
        end
        IPPVDm = TEC2m(IPPVD);
        a = unique(IPPVDm(:,1));
        startdata = find(IPPVDm(:,1) == a(1), 1, 'first');
        enddata = size(IPPVDm,1);

        SPF  = PFmodel12_fast(IPPVDm, 1, 3, [startdata, enddata]);
        
        prepDir = fullfile(workspace, 'preprocessing');
        if exist(prepDir, 'dir') == 0
            mkdir(prepDir);
        end
        saveFile = fullfile(prepDir, 'PFresult.mat');
        parsave(saveFile, 'SPF', SPF);
        ok(idx) = true;

    catch ME
        ok(idx) = false;
        st = ME.stack;
        if ~isempty(st)
            top = sprintf("%s (line %d)", st(1).name, st(1).line);
        else
            top = "no stack";
        end
        errLog{idx} = struct( ...
            "k", k, ...
            "workspace", workspace, ...
            "identifier", ME.identifier, ...
            "message", ME.message, ...
            "top", top, ...
            "stack", st );
    end
end
failed = find(~ok);
fprintf("\nDone. Success=%d, Failed=%d\n", sum(ok), numel(failed));
for ii = 1:numel(failed)
    e = errLog{failed(ii)};
    fprintf("\n[FAIL] k=%d\nworkspace: %s\n%s\n%s: %s\n", ...
        e.k, e.workspace, e.top, e.identifier, e.message);
end
save(fullfile(base_path, 'parfor_error_log.mat'), 'ok', 'wsLog', 'errLog');



parfor k =19:30
    % 设置当前文件夹的 workspace 路径
    workspace = fullfile(base_path, folders(k).name);
    ippcomplete=dir(fullfile(workspace,'e_s_vD_ipp'));
    if numel(ippcomplete)==3
        A=load(fullfile(workspace,'e_s_vD_ipp','vD_all.mat'));
        IPPVD=A.IPPVD;
    else
    IPPVD=extractIPPVD(workspace);%unit:TECU
    end
    % IPPVD=IPPVD_data{k-7};
    IPPVDm=TEC2m(IPPVD); % into unit:m，经度范围-180到180
    a=unique(IPPVDm(:,1)); %已经按照时间排序
    startdata = find(IPPVDm(:,1) == a(1), 1, 'first');
    enddata=size(IPPVDm,1);
    SPF=PFmodel12_fast(IPPVDm, 1, 3,[startdata,enddata]); %仅去除空间趋势项，方便做S+S，S+ST
    STPF=PFmodel12_fast(IPPVDm, 2, 3,[startdata,enddata]);%去除时间与空间趋势项，做ST+ST
    if exist([workspace '\preprocessing'],'dir')==0
        mkdir([workspace '\preprocessing']);
    end
    parsave([workspace '\preprocessing\PFresult.mat'],'SPF',SPF,'STPF',STPF);
end

%% Step2: Collect the daily spatial exvariograms and temporal exvariograms
%由于实时，因此不需要worksapce存储数据，直接导出离散值即可，按照epoch迭代variogram
% Model 1：单S模型，被抛弃
% Grids definition
% [LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
% IGP_list=[LatC(:), LonC(:)];
% 
% base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\Hongkong\download';
% folders = dir(base_path);
% folders = folders([folders.isdir]); % 只保留目录
% folders = folders(~ismember({folders.name}, {'.', '..'})); % 排除 '.' 和 '..'
% for k = 1:length(folders)
% workspace = fullfile(base_path, folders(k).name);
% load([workspace '\preprocessing\PFresult.mat']);
% S_S=variogram_spatial([SPF(:,3) SPF(:,4)  SPF(:,1)] ,SPF(:,7),'sth',400);% dist,vario,epoch
% S_T=variogram_temporal([SPF(:,3) SPF(:,4)  SPF(:,1) SPF(:,7)] ,'sth',2);% time diff,vario,pre epoch,now epoch
% % Model 2: ST模型得到的variogram
% ST_S=variogram_spatial([STPF(:,3) STPF(:,4)  STPF(:,1)] ,STPF(:,7),'sth',400);
% ST_T=variogram_temporal([STPF(:,3) STPF(:,4)  STPF(:,1) STPF(:,7)],'sth',2);
% % S和ST模型下的离散ST semvariogarms
% S_ST=variogram_ST([SPF(:,3) SPF(:,4)  SPF(:,1)] ,SPF(:,7),'sth',400); %数据量太大一个历元香港地区有30000个点对是否需要bin呢
% ST_ST=variogram_ST([STPF(:,3) STPF(:,4)  STPF(:,1)] ,STPF(:,7),'sth',400);
%     if exist([workspace '\preprocessing'],'dir')==0
%         mkdir([workspace '\preprocessing']);
%     end
%     save([workspace '\preprocessing\EX_vario.mat'],'S_S','S_T','ST_S','ST_T','S_ST','ST_ST');
% end

%% Step 3: Construct Real time ST variogram models at IGPs
params = struct('sth_km',2100,'win_deg',15,'time_back',300,'time_step',30,...
    'use_single',true,'code_version','v1.2');
[LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
IGP_list=[LatC(:), LonC(:)];
for k =14:30
    outdir = fullfile(base_path, folders(k).name, 'preprocessing');
    outdir3 = fullfile(base_path, folders(k).name, 'NLSex');
    outdir2 = fullfile(base_path, folders(k).name, 'NLSprocessing');

    if exist(outdir3,'dir')==0
        mkdir(outdir3);
    end
    if exist(outdir2,'dir')==0
        mkdir(outdir2);
    end
    S = load(fullfile(outdir, 'PFresult.mat'), 'STPF');
    % 减少广播：每个 worker 一次性获取
    STPF_const = parallel.pool.Constant(S.STPF);
  
    parfor i = 1:78
        t = getCurrentTask(); wid = 0; if ~isempty(t), wid = t.ID; end
        if mod(i,10)==1 
        fprintf('Worker %d handling k=%d at IGP=%d\n', wid, k, i);
        end
        STPF = STPF_const.Value;  % 取共享数据
        save_IGP_variogram_fast3(STPF, IGP_list(i,1), IGP_list(i,2), params, outdir3);
    end

     parfor i = 1:78
        t = getCurrentTask(); wid = 0; if ~isempty(t), wid = t.ID; end
        if mod(i,10)==1 
        fprintf('Worker %d handling k=%d at IGP=%d\n', wid, k, i);
        end
        STPF = STPF_const.Value;  % 取共享数据
        save_IGP_variogram_fast2(STPF, IGP_list(i,1), IGP_list(i,2), params, outdir2);
     end

      parfor i = 1:78
        t = getCurrentTask(); wid = 0; if ~isempty(t), wid = t.ID; end
        if mod(i,10)==1 
        fprintf('Worker %d handling k=%d at IGP=%d\n', wid, k, i);
        end
        STPF = STPF_const.Value;  % 取共享数据
        save_IGP_variogram_fast(STPF, IGP_list(i,1), IGP_list(i,2), params, outdir);
    end
    % （可选）释放 Constant 以便下一轮外层 for 的数据能被回收
    delete(STPF_const);
end
%% Step 4: Generate and GIVD testing on each IGP
base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\USA Quiet2';
% 获取所有文件夹列表
folders = dir(base_path);
folders = folders([folders.isdir]); % 只保留目录
folders = folders(~ismember({folders.name}, {'.', '..'})); % 排除 '.' 和 '..'
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('进度更新: %s\n', msg));
[LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
IGP_list=[LatC(:), LonC(:)];
for k=14:15
    %load MT18：GIVD and GIVE
    workspace = fullfile(base_path, folders(k).name);
     outdir = fullfile(base_path, folders(k).name, 'preprocessing');
     outdir2=fullfile(base_path, folders(k).name, 'NLSprocessing');
     outdir3=fullfile(base_path, folders(k).name, 'NLSex');
    givd=load(fullfile(workspace,'Broadcast\Broadcast GIVD.mat'));%0s到86370s
    givd=givd.givd;
    % load('E:\PolyU_PhD\spatialvario\21308\Broadcast\Broadcast GIVE.mat');
    igpmask_tmp=load(fullfile(workspace,'Broadcast\Broadcast IGPs.mat'));
    igpmask_tmp=igpmask_tmp.igpmask_tmp;

    mat_file = fullfile(workspace, 'preprocessing', 'PFresult.mat');
    data = load(mat_file); % 加载 PFresult.mat

    if exist([fullfile(workspace) '\NLSGIVD4'],'dir')==0
        mkdir([fullfile(workspace) '\NLSGIVD4']);
    end
    if exist([fullfile(workspace) '\Simulate GIVD'],'dir')==0
        mkdir([fullfile(workspace) '\Simulate GIVD']);
    end
    %load IGP-i variogram model
    % 初始化并行池（可选，视硬件调整线程数）
    % 创建 DataQueue 用于进度监控
    parfor i=1:size(IGP_list,1)
        PSOMODELs=load([outdir '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel=PSOMODELs.ST_STmodels;
        NLSMODELs=load([outdir2 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel2=NLSMODELs.ST_STmodels;
        ExMODELs=load([outdir3 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel3=ExMODELs.ST_STmodels;

        [A,~] = NLSSTGIVx2([800,2100], [10,40], data.STPF, 1, ST_STmodel,ST_STmodel2, ST_STmodel3,'GIVD',IGP_list(i,:), givd, igpmask_tmp);
        outfile = fullfile(workspace,'NLSGIVD4', sprintf('IGP_%s_%s_STKGIVD.mat', num2str(IGP_list(i,1)), num2str(IGP_list(i,2))));
       parsave(outfile, 'A',A); 
        [C,~] = UnKGIVx([800,2100], [10,30], data.SPF, 1, 'GIVD',IGP_list(i,:), givd, igpmask_tmp); 
        outfile2 = fullfile(workspace,'Simulate GIVD', sprintf('IGP_%s_%s_UnKGIVD.mat', num2str(IGP_list(i,1)), num2str(IGP_list(i,2))));
        parsave(outfile2, 'C',C);
    end
    %分析部分转入GIVD compare.mat
    % [A,B] = STGIVx2([800,2100], [10,30], PF_data, 1, ST_STmodels, 'GIVD',IGPs, givd, igpmask_tmp);
    % STGIVx2 为用UnK补偿的方案
    % [C,D] = UnKGIVx([800,2100], [10,30], PF_data, 1, 'GIVD',IGPs, givd, igpmask_tmp);
end
%% Step 5:Generate the raw threat model using STK and UnK
for k=1:13
    %load MT18：GIVD and GIVE
    workspace = fullfile(base_path, folders(k).name);
    outdir = fullfile(base_path, folders(k).name, 'preprocessing');
    givd=load(fullfile(workspace,'Broadcast\Broadcast GIVD.mat'));%0s到86370s
    givd=givd.givd;
    % load('E:\PolyU_PhD\spatialvario\21308\Broadcast\Broadcast GIVE.mat');
    igpmask_tmp=load(fullfile(workspace,'Broadcast\Broadcast IGPs.mat'));
    igpmask_tmp=igpmask_tmp.igpmask_tmp;
    mat_file = fullfile(workspace, 'preprocessing', 'PFresult.mat');
    data = load(mat_file); % 加载 PFresult.mat
    STthreat_parts = cell(size(IGP_list,1),1);
    parfor i=1:size(IGP_list,1)
        MODELs= load([outdir '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel=MODELs.ST_STmodels;
        % % IGPthreat=STthreat2(ST_STmodels,[800,2100], [10,30], data.STPF, 1,6,IGP_list(i,:), givd, igpmask_tmp);
        IGPthreat=STIPPthreat_fast(ST_STmodel,[800,2100], [10,30], data.STPF, 1,5,IGP_list(i,:), givd, igpmask_tmp,dq);
        STthreat_parts{i} = IGPthreat; 
        % IGPUnK=simplethreat2([800,2100], [10,30], data.SPF, 1, 200,IGP_list(i,:), givd, igpmask_tmp);
        % UnKIGPthreat=[UnKIGPthreat;IGPUnK];
    end
    STthreat2 = vertcat(STthreat_parts{:}); 
    STthreat=STthreat2(STthreat2(:,16)>0,:);
    UnKIPPthreat=simplethreat_fast([800,2100], [10,30], data.SPF, 1, 200,dq);
    UnKthreat=UnKIPPthreat(UnKIPPthreat(:,18)>0,:);
    if exist([workspace '\threat model construction'],'dir')==0
        mkdir([workspace '\threat model construction']);
    end
    save(fullfile(workspace,'threat model construction\STthreat_all.mat'),'STthreat2','-v7.3');
    save(fullfile(workspace,'threat model construction\UnKthreat_all.mat'),"UnKIPPthreat",'-v7.3');
    save(fullfile(workspace,'threat model construction\STthreat_valid.mat'),"STthreat",'-v7.3');
    save(fullfile(workspace,'threat model construction\UnKthreat_valid.mat'),"UnKthreat",'-v7.3');
end

%% Generate the raw threat model using STK and UnK,但是统一由IPP生成
base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\USA';
% 获取所有文件夹列表
folders = dir(base_path);
folders = folders([folders.isdir]); % 只保留目录
folders = folders(~ismember({folders.name}, {'.', '..'})); % 排除 '.' 和 '..'
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('进度更新: %s\n', msg));
[LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
IGP_list=[LatC(:), LonC(:)]; 
for k=1:13
    workspace = fullfile(base_path, folders(k).name);
  outdir = fullfile(base_path, folders(k).name, 'preprocessing');
     outdir2=fullfile(base_path, folders(k).name, 'NLSprocessing');
     outdir3=fullfile(base_path, folders(k).name, 'NLSex');


    mat_file = fullfile(workspace, 'preprocessing', 'PFresult.mat');
    data = load(mat_file); % 加载 PFresult.mat
    STthreat_parts = cell(size(IGP_list,1),1);
    parfor i=1:size(IGP_list,1)
        PSOMODELs=load([outdir '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel=PSOMODELs.ST_STmodels;
        NLSMODELs=load([outdir2 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel2=NLSMODELs.ST_STmodels;
        ExMODELs=load([outdir3 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel3=ExMODELs.ST_STmodels;

        IGPthreat=STIPPthreat_fast_IPPbox( ST_STmodel,ST_STmodel2, ST_STmodel3,[800,2100], [10,40], data.STPF, 1,5,IGP_list(i,:), dq);
        STthreat_parts{i} = IGPthreat; 
    end
    STthreat2 = vertcat(STthreat_parts{:}); 
    STthreat=STthreat2(STthreat2(:,16)>0,:);
    if exist([workspace '\threat model construction2'],'dir')==0
        mkdir([workspace '\threat model construction2']);
    end
    parsave(fullfile(workspace,'threat model construction2\STthreat_all.mat'),'STthreat2',STthreat2);
    parsave(fullfile(workspace,'threat model construction2\STthreat_valid.mat'),"STthreat",STthreat);
end


%转入constructthreat进行分析
%% Step 6: Fusion threat and construct threat model
% Go To：4Dplot.m
STthreat_parts=cell(numel(folders),1);
UnKthreat_parts=cell(numel(folders),1);
for k=1:13
    %load MT18：GIVD and GIVE
    workspace = fullfile(base_path, folders(k).name);
    mat_file = fullfile(workspace, 'threat model construction2', 'STthreat_valid.mat');
    mat_file2 = fullfile(workspace, 'threat model construction', 'UnKthreat_valid.mat');
    data = load(mat_file); % 加载 PFresult.mat
    STthreat_parts{k} = data.STthreat;
    data2= load(mat_file2);
    UnKthreat_parts{k}=data2.UnKthreat;
end
STthreat_validfusion=vertcat(STthreat_parts{:});
UnKthreat_validfusion=vertcat(UnKthreat_parts{:});


%% GIVE Analysis-
% Obj1: GIVE reduction at IGPs all results;
% Obj2: GIVE<=6m percentage,each IGP at竖状图，and reduction
% Obj3: GIVE and error bin counts
% Obj4: sigma_undersampled/sigma_GIVE reduction
A_give=givd2give(A,12,13,14,15,11,sigma_mono_4D);
C_give=givd2give(C,12,13,11,sigma_mono);

result=[mean(A_give(:,18:19)),mean(C_give(:,18:19))];
mean(A_give(:,18)./A_give(:,19))
mean(C_give(:,18)./C_give(:,19))
max(A_give(:,19))
max(C_give(:,19))

%% UIVD and UIVE
%Firstly,we should estimate the GIVD values from STK and UnK:坏处是有可能GIVD与观测量不一致 
% Or we estimate the user directly without
% GIVDs：坏处是STK用哪里的variogram去拟合，用IGP则需要每次都引入
base_path = 'I:\主机材料_polyu PhD\Spatial Temporal Threat Model\USA';
% 获取所有文件夹列表
folders = dir(base_path);
folders = folders([folders.isdir]); % 只保留目录
folders = folders(~ismember({folders.name}, {'.', '..'})); % 排除 '.' 和 '..'
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('进度更新: %s\n', msg));
[LatC, LonC]=meshgrid(25:5:50,-130:5:-70);
IGP_list=[LatC(:), LonC(:)]; 
for k=1:13
    workspace = fullfile(base_path, folders(k).name);
    outdir = fullfile(base_path, folders(k).name, 'preprocessing');
    outdir2=fullfile(base_path, folders(k).name, 'NLSprocessing');
    outdir3=fullfile(base_path, folders(k).name, 'NLSex');
    mat_file = fullfile(workspace, 'preprocessing', 'PFresult.mat');
    data = load(mat_file); % 加载 PFresult.mat
    STUIVD_parts = cell(size(IGP_list,1),1);
    %load IGP-i variogram model
    % 初始化并行池（可选，视硬件调整线程数）
    % 创建 DataQueue 用于进度监控
    parfor i=1:size(IGP_list,1)
        PSOMODELs=load([outdir '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel=PSOMODELs.ST_STmodels;
        NLSMODELs=load([outdir2 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel2=NLSMODELs.ST_STmodels;
        ExMODELs=load([outdir3 '\IGP_' num2str(IGP_list(i,1)) 'N_' num2str(abs(IGP_list(i,2))) 'W_variogram.mat']);
        ST_STmodel3=ExMODELs.ST_STmodels;
        STUIVD_parts{i} = STIPP_fast(ST_STmodel,ST_STmodel2, ST_STmodel3,[800,2100], [10,40], data.STPF, 1,5,IGP_list(i,:), dq);
    end
    STKUIVD= vertcat(STUIVD_parts{:});
      if exist([workspace '\UIVD'],'dir')==0
        mkdir([workspace '\UIVD']);
    end
    parsave(fullfile(workspace,'UIVD\STKUIVD.mat'),'STKUIVD',STKUIVD);
end

for k=1:13
    workspace = fullfile(base_path, folders(k).name);
    mat_file = fullfile(workspace, 'preprocessing', 'PFresult.mat');
    data = load(mat_file); % 加载 PFresult.mat
    UnKUIVD=UnKIPP_fast([800,2100], [10,30], data.SPF, 1, dq);
    parsave(fullfile(workspace,'UIVD\UnUIVD.mat'),'UnUIVD',UnKUIVD);
end
 % outfile2 = fullfile(workspace,'Simulate GIVD', sprintf('IGP_%s_%s_UnKGIVD.mat', num2str(IGP_list(i,1)), num2str(IGP_list(i,2))));
 % parsave(outfile2, 'C',C);

