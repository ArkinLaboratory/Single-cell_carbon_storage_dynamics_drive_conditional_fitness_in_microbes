clear;

% folder='F:\both_inc_dec_2_to_0dot1_new_cropped\inc\summary_sim';

% folder='F:\Experiment-346_change_ion\lane1\summary_sim';

folder='F:\change_ion_third\lane2\summary_sim';


% folder='E:\both_inc_dec_2_to_0dot1_new_cropped\inc\summary_sim';

% folder='E:\Experiment-299_Lane1_01_001_0001_Lane2_1_01_001-V7\lane1\summary_sim';

postprocfolder=[folder,filesep];

% start_time=1;
% end_time=900;
% smedia=[0,262,476,697];

%change 4 con
start_time=1;
end_time=997;
% end_time=340;



start_time=1;
end_time=997;

thres_bright = 300;
thres_dark = 100;
areathres = 400;
areathresup = 800;
minsize=0;

tic;
drawind = 0;
pp=1;
thres_gen = 50;

sample_pixels=cell(end_time-start_time+1,1);
sample_granules=cell(end_time-start_time+1,1);
sample_granules_pixel=cell(end_time-start_time+1,1);
sample_area=cell(end_time-start_time+1,1);
sample_rate_smooth=cell(end_time-start_time+1,1);
sample_rate=cell(end_time-start_time+1,1);
sample_divs=zeros(end_time-start_time+1,1);

sample_outside_gran_without=cell(end_time-start_time+1,1);
sample_outside_gran=cell(end_time-start_time+1,1);
sample_totf=cell(end_time-start_time+1,1);
sample_moving_outside_gran=cell(end_time-start_time+1,1);
sample_moving_totf=cell(end_time-start_time+1,1);
sample_traj_high=cell(end_time-start_time+1,1);


sample_totf_outside_gran=cell(end_time-start_time+1,1);
sample_dilute=cell(end_time-start_time+1,1);
max_gran_size=10000;
rate_granule_size=cell(1,max_gran_size);
smooth_rate_granule_size=cell(1,max_gran_size);

darknum=0;
darknum1=0;
outn=0;
res1=[];

stat_cell_highf_frames = [];


dd=dir(postprocfolder);
output=[];
chamber = 0;

moving_ave=20;



draw_traj=0;
traj_len=100;



for idir=1:length(dd)
    if length(dd(idir).name)>6 && strcmp(dd(idir).name(1:6),'Trench')
        load([postprocfolder,'\',dd(idir).name])
        
        chamber = chamber+1;
        
        display(dd(idir).name);
        for val=1:length(lineage)
            
            ll=0;
            allframe=length(lineage{val}.framenbs);
            framenbs = lineage{val}.framenbs;
            st=find(framenbs>=start_time,1);
            et=find(framenbs>end_time,1);
            
            if isempty(et)
                et=allframe;
            else
                et=et-1;
            end
            
            fluo = lineage{val}.fluo1;
            if et>length(lineage{val}.length) || et>length(fluo)
                et=min(length(fluo),length(lineage{val}.length));
            end
            
            if ~isempty(st) && et>=st
                
                
                localFramenbs=framenbs(st:et);
                
                y=lineage{val}.length(st:et);
                
                ia=0;
                area=zeros(et-st+1,1);
                for i=st:et
                    ia=ia+1;
                    if length(size(lineage{val}.pixel_xy)) == 3
                        area(ia) = size(lineage{val}.pixel_xy,3);
                    else
                        area(ia) = length(lineage{val}.pixel_xy{i,1});
                    end
                end
                
                ia=0;
                gran=zeros(et-st+1,1);
                mgran=zeros(et-st+1,1);
                pixels=cell(et-st+1,1);

                outgran=zeros(et-st+1,1);
              
                for i=st:et
                    ia=ia+1;
                    if iscell(fluo)
                        gran(ia) = sum(fluo{i}>thres_bright);
                        mgran(ia) = max(fluo{i});
                        outgran(ia) = mean(fluo{i}(fluo{i}<=thres_bright));
                    else
                        gran(ia) = sum(fluo(i,:)>thres_bright);
                        mgran(ia) =  max(fluo(i,:));
                        curr_fluo = fluo(i,:);
                        outgran(ia) = mean(curr_fluo(curr_fluo<=thres_bright));
                    end
                end
                
                
                dau = lineage{val}.daughters(st:et);
                
%                 divs=find(dau);
%                 sf=1;
%                 
%                 for i=1:length(divs)
%                     ef=divs(i);
%                     
%                     ff=localFramenbs(ef)-start_time+1;
%                     sample_divs(ff) = sample_divs(ff)+1;
%                 end
                
                for i=1:length(gran)
                    if mgran(i)>=thres_dark && area(i) > areathres
                        ff=localFramenbs(i)-start_time+1;
                        
                        out=[sample_granules{ff},gran(i)/area(i)];
                        sample_granules{ff}=out;
                        

                        out=[sample_area{ff},area(i)];
                        sample_area{ff}=out;

                        out=[sample_outside_gran{ff}, outgran(i)];
                        sample_outside_gran{ff} = out;
                        
                        if dau(i)~=0
                            sample_divs(ff) = sample_divs(ff)+1;
                        end
                        
                        if mgran(i) < thres_bright
                            sample_outside_gran_without{ff} = sample_outside_gran{ff};
                        end

                    else
                        darknum=darknum+1;
                    end
                end

                if drawind && outn<50 && et-st+1 > 500
                    outn=outn+1;
                    figure('Renderer', 'painters', 'Position', [10 20 2000 300])
                    hold on;
                    yyaxis left
                    
                    for t=2:length(y)
                        if abs(y(t)-y(t-1))/y(t-1)>0.4 && t<length(y) && dau(t+1)>0
                            y(t)=y(t-1);
                            gran(t)=gran(t-1)
                            %y(t)=y(t-1)+abs(y(t-1)-y(t-2));
                            if t<length(y)
                                dau(t+1)=0;
                            end
                        end
                        
                    end
                    
                    plot(localFramenbs,y*0.055,'-',localFramenbs(dau>0), y(dau>0)*0.055, 'ro')
                    
                    xlim([localFramenbs(1) localFramenbs(end)])
                    ylim([0 6])
                    ylabel('Cell Length (\mum)')
                    
                    yyaxis right
                    %
                    
                    %plot(localFramenbs,gran,'b-',localFramenbs(1:length(outgran)),outgran,'r-');
                    plot(localFramenbs,gran,'-');
                    ylim([-10 1000])
                    ylabel('# bright pixels')
                    xticks(localFramenbs(1)-1:40:localFramenbs(end))
                    xlabel('Time (frames per 5mins)')

                    set(gca,'FontName','Times new roman','FontSize',15)
                    title(['pos=',num2str(1),' trench=',num2str(chamber),' val=',num2str(val),' # div=',num2str(sum(dau>0))])
                    box on;
                end
            end
            
        end
    end
end
display(darknum);
toc;


%%
stat=zeros(length(sample_granules),2);
numcell=zeros(length(sample_granules),1);
for i=1:length(sample_granules)
    stat(i,1)=mean(sample_granules{i});
    stat(i,2)=var(sample_granules{i});
    numcell(i)=length(sample_granules{i});
end

statarea=zeros(length(sample_area),2);
for i=1:length(sample_area)
    statarea(i,1)=mean(sample_area{i});
    statarea(i,2)=var(sample_area{i});
end

mean_outgran = zeros(length(sample_outside_gran), 2);
for i = 1:length(sample_outside_gran)
    raw_outgran = sample_outside_gran{i};
    no_nan_outgran = raw_outgran(~isnan(raw_outgran));
    mean_outgran(i, 1) = mean(no_nan_outgran);
    mean_outgran(i, 2) = var(sample_outside_gran{i});
end

mean_without = zeros(length(sample_outside_gran_without), 2);
for i = 1:length(sample_outside_gran_without)
    raw_without = sample_outside_gran_without{i};
    no_nan_without = raw_without(~isnan(raw_without));
    mean_without(i, 1) = mean(no_nan_without);
    mean_without(i, 2) = var(sample_outside_gran_without{i});
end
%% 
% figure;
% hold on
% yyaxis left
% 
% plot(start_time:end_time,stat(:,1),'-','linewidth',1)
% 
% set(gca,'FontName','Times new roman','FontSize',15,'linewidth',1)
% ylim([0,0.20])
% ylabel('Fraction')
% 
% xlim([start_time,end_time]);
% xlabel('Time (frame, 10min)')
% 
% yyaxis right
% 
% set(gca,'FontName','Times new roman','FontSize',15)
% plot(start_time:end_time, mean_outgran(:,1), '-g', 'linewidth', 1);
% plot(start_time:end_time, mean_without(:,1),'--','Color', [0, 0.5, 0], 'linewidth', 1);
% ylim([120,260])
% 
% % con = [0.1,0.5,1,2];
% % con=[2,1,0.5,0.1];
% smedia=[0,205,415,636];
% con = {'0.1', '0.1', '0.1',  '0.1'};
% for i=1:length(smedia)
%     xline(smedia(i),'linewidth',1)
%     txt = con{i};
%     text(smedia(i),250,txt,'HorizontalAlignment','left')
% end
% 
% legend('Mean PHB granule over cell area', 'outgran', 'without');
% 
% yyaxis left
% 
% set(gca,'FontName','Times new roman','FontSize',15)

%%
figure;
hold on
yyaxis left

plot(start_time:end_time,stat(:,1),'r-')

set(gca,'FontName','Times new roman','FontSize',15,'linewidth',1)
ylim([0,0.20])
ylabel('Fraction')
% set(gca,'yticklabel',{'0%', '5%', '10%', '15%' ,'20%'});
plot(start_time:end_time,sample_divs./numcell,'Color',[0.7, 0.7, 0])


xlim([start_time,end_time]);
xlabel('Time (frame, 10min)')

yyaxis right

% plot(start_time:end_time,numcell,'linewidth',1)

set(gca,'FontName','Times new roman','FontSize',15)
plot(start_time:end_time,statarea(:,1),'Color',[0, 0.5, 0])

ylim([300,1300])
smedia=[0,205,415,636];
con = {'0.1N', '0.5N', '1N',  '2N'};


for i=1:length(smedia)
    xline(smedia(i),'linewidth',1)
    txt = con{i};
    text(smedia(i),1270,txt,'HorizontalAlignment','left')
end

yyaxis left
% plot([start_time,end_time], [0.052, 0.052], '-.')
% plot([start_time,end_time], [0.065, 0.065], '-.')
legend('Mean PHB granule over cell area', 'Division events per cell','Cell Area')
title('ion constant change N')
set(gca,'FontName','Times new roman','FontSize',15)


%%
% clear res;
% save(['change_c_',num2str(thres_bright),'.mat']);