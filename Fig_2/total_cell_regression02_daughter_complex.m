clear;
start_time=1;
end_time=961;

thres_bright = 300;
thres_dark = 80;
areathres = 400;
areathresup = 800;
minsize=0;

% folder='D:\both_inc_dec_2_to_0dot1_new\dec\summary_sim';
% folder='D:\both_inc_dec_2_to_0dot1_new\inc\summary_sim';
% postprocfolder=[folder,filesep];
% load mothers4-20.mat;
% load mothers51-60.mat;

%% 

% parpool(16)


%%
needplot=0;
numlines=0;
thres=80;
thres_r2 = 0.8;
ave=0;
totrate=[];
totout=[];

rates=[];

cell_data_collect = {};
total_daughter01_data_collect = [];
total_daughter02_data_collect = [];
% time_cell02 = [];

% for num_mother = 1 : 12
% for num_mother = 25

% dd=dir(postprocfolder);
output=[];
chamber = 0;

% for fov=[4]
% for fov=[4:20]


for fov=[41:60]
% for fov=41:43

    for crop=1:3
    % for crop=2
        name=['pos',num2str(fov),'_',num2str(crop)];

%         postprocfolder = ['F:\both_inc_dec_2_to_0dot1_new_cropped\inc\output1_',name];
%         postprocfolder = ['F:\both_inc_dec_2_to_0dot1_new_cropped\dec\output1_',name];
        postprocfolder = ['F:\both_inc_dec_2_to_0dot1_new_cropped\dec\output1_',name];

        dd=dir(postprocfolder);
%         for idir=3:4 
        for idir=1:length(dd)

            if length(dd(idir).name)>6 && strcmp(dd(idir).name(1:6),'Trench')
                load([postprocfolder,'\',dd(idir).name])
                lineage = res.lineage;

                
                for num_mother = 1 : length(lineage)
%                 for num_mother = 1
                    num_nonzero = nnz(lineage{1, num_mother}.daughters);
        
                    % 如果非零值数量 < 2，则跳过
                    if num_nonzero < 2  
                        continue; % 直接跳过该 i，进入下一个 i
                    end
                    % disp(['FOV: ', num2str(fov), num2str(crop)]);


                    % cmother=lineage(1, num_mother);
                    % cmother = cmother{1, 1};
                    cmother=lineage{1, num_mother};

                    this_logic = logical(cmother.daughters(:));
                
                    divisions = find(this_logic);

                    for t = 1:length(divisions)
                        
                        if divisions(t) == 1 || t+1 > length(divisions)
                            break;
                        end
                        
                        if t == 1
                            prev_time = 1;
                        else
                            prev_time = divisions(t-1);
                        end
                        log_size = log(cmother.length(1,prev_time:divisions(t)-1));
                        
                        if length(log_size) < 6
                            continue
                        end
                        linear_regres = fitlm(prev_time:divisions(t)-1, log_size);
                        r2 = linear_regres.Rsquared.Ordinary;
                        
                        if isnan(r2)
                            continue
                        end
                        if r2 < thres_r2
                            continue
                        end

                        % daughter 1
                        prev_time = divisions(t);

                        log_size = log(cmother.length(1,prev_time:divisions(t+1)-1));
                        
                        if length(log_size) < 6
                            continue
                        end
                        linear_regres = fitlm(prev_time:divisions(t+1)-1, log_size);
                        r2 = linear_regres.Rsquared.Ordinary;
                        
                        if isnan(r2)
                            continue
                        end
                        if r2 < thres_r2
                            continue
                        end

                        % daughter 2
                        num_daughter = cmother.daughters(divisions(t));
                        cdaughter = lineage{1, num_daughter};

                        % 如果 daughters 里面没有至少一个非零值，跳过该 i
                        if isempty(cdaughter.daughters) || all(cdaughter.daughters == 0) 
                            continue; % 跳过当前 i，进入下一个 i
                        end
                        this_logic_dau = logical(cdaughter.daughters(:));
                
                        divisions_dau = find(this_logic_dau);
                        div_index_dau = 1;
                               
                        if divisions_dau(1) == 1
                            break;
                        end                       

                        log_size = log(cdaughter.length(1,1:divisions_dau(1)-1));
                        
                        if size(log_size) < 6
                            continue
                        end
                        linear_regres = fitlm(1:divisions_dau(1)-1, log_size);
                        r2 = linear_regres.Rsquared.Ordinary;
                        
                        if isnan(r2)
                            continue
                        end
                        if r2 < thres_r2
                            continue
                        end


                        % daughter 1 PHB
                        gran=[];
                        fluo_raw=cmother.fluo1;
                        for j=divisions(t):divisions(t+1)
            
                            if iscell(fluo_raw)
                                try 
                                    fluo = fluo_raw{j};
                                catch 
                                    continue;
                                end
                            else
                                try 
                                    fluo = fluo_raw(j,:);
                                catch
                                    continue;
                                end
                            end
                            gran(j - prev_time + 1) = sum(fluo>thres_bright);
                            % gran(j) = sum(fluo>thres_bright);
                        end
                        display("dau1", num2str(divisions(t)))
                        mdl = fitlm(1:length(gran), gran(:));
                        prod_rate_dau1 = mdl.Coefficients.Estimate(2);

                        pos_dau1 = mean(cmother.cellpixels{divisions(t),2});

                        % mother PHB
            
                        if iscell(fluo_raw)
                            try 
                                fluo = fluo_raw{divisions(t)-1};
                            catch 
                                continue;
                            end
                        else
                            try 
                                fluo = fluo_raw(divisions(t)-1,:);
                            catch
                                continue;
                            end
                        end
                        
                        mother_gran = sum(fluo > thres_bright);
                        display("mother", num2str(mother_gran))

                        % daughter 2 PHB

                        gran_dau=[];
            
                        for j=1:divisions(1)
            
                            fluo=cdaughter.fluo1;
                            if iscell(fluo)
                                try 
                                    fluo = fluo{j};
                                catch 
                                    continue;
                                end
                            else
                                try 
                                    fluo = fluo(j,:);
                                catch
                                    continue;
                                end
                            end
                            gran_dau(j) = sum(fluo>thres_bright);
                        end
                        display("dau2", num2str(length(gran_dau)))
                        mdl = fitlm(1:length(gran_dau), gran_dau(:));
                        prod_rate_dau2 = mdl.Coefficients.Estimate(2);

                        pos_dau2 = mean(cdaughter.cellpixels{1,2});

                        end_data1 = size(total_daughter01_data_collect, 2);
                        total_daughter01_data_collect(1, end_data1+1) = size(cmother.fluo1{1,divisions(t)},2);
                        total_daughter01_data_collect(2, end_data1+1) = size(cmother.fluo1{1,divisions(t+1)-1},2);
                        total_daughter01_data_collect(3, end_data1+1) = gran(1);
                        total_daughter01_data_collect(4, end_data1+1) = gran(end);
                        total_daughter01_data_collect(5, end_data1+1) = gran(end)-gran(1);
                        total_daughter01_data_collect(6, end_data1+1) = gran(1)/total_daughter01_data_collect(2, end_data1+1);
                        total_daughter01_data_collect(7, end_data1+1) = 1/(divisions(t+1)-divisions(t));
                        total_daughter01_data_collect(8, end_data1+1) = divisions(t+1)-divisions(t);
                        total_daughter01_data_collect(9, end_data1+1) = prod_rate_dau1;
                        total_daughter01_data_collect(10, end_data1+1) = gran(1)/mother_gran;
                        total_daughter01_data_collect(11, end_data1+1) = pos_dau1;

                        end_data2 = size(total_daughter02_data_collect, 2);
                        total_daughter02_data_collect(1, end_data2+1) = size(cdaughter.fluo1{1,1},2);
                        total_daughter02_data_collect(2, end_data2+1) = size(cdaughter.fluo1{1,divisions_dau(1)-1},2);
                        total_daughter02_data_collect(3, end_data2+1) = gran_dau(1);
                        total_daughter02_data_collect(4, end_data2+1) = gran_dau(end);
                        total_daughter02_data_collect(5, end_data2+1) = gran_dau(end)-gran_dau(1);
                        total_daughter02_data_collect(6, end_data2+1) = gran_dau(1)/total_daughter02_data_collect(2, end_data2+1);
                        total_daughter02_data_collect(7, end_data2+1) = 1/(divisions_dau(1));
                        total_daughter02_data_collect(8, end_data2+1) = divisions_dau(1);
                        total_daughter02_data_collect(9, end_data2+1) = prod_rate_dau2;
                        total_daughter02_data_collect(10, end_data2+1) = gran_dau(1)/mother_gran;
                        total_daughter02_data_collect(11, end_data2+1) = pos_dau2;
                        
    

                    end
                    
                end
            end
        end
    end
end

%% 
% 计算数据差值
hist_data = total_daughter02_data_collect(11,:) - total_daughter01_data_collect(11,:);

% 创建直方图，设置步长为 0.5
figure;
histogram(hist_data, 'BinWidth', 0.5, 'FaceColor', [0.5 0.5 0.5]); % 灰色柱子

% 添加标题和标签
title('Histogram of Differences');
xlabel('Difference in (11,:) Values');
ylabel('Frequency');

% 显示网格以便更清晰
grid on;

total_daughter02_data_collect(12,:) = hist_data;


% 提取数据
data_x = total_daughter02_data_collect(10,:);
grouping_data = total_daughter02_data_collect(12,:);

% 按照 total_daughter02_data_collect(12,:) 进行分类
data_positive = data_x(grouping_data > 0); % 对应正数的
data_negative = data_x(grouping_data < 0); % 对应负数的

% 创建直方图
figure;
hold on; % 保持同一张图

% 绘制正数组数据的直方图（蓝色）
histogram(data_positive, 'BinWidth', 0.1, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Positive');

% 绘制负数组数据的直方图（红色）
histogram(data_negative, 'BinWidth', 0.1, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Negative');

% 添加标题和轴标签
title('Histogram of total\_daughter02\_data\_collect(10,:)');
xlim([0,1.5])
xlabel('Value');
ylabel('Frequency');

% 添加图例
legend('show');

% 显示网格
grid on;

% 关闭 hold
hold off;