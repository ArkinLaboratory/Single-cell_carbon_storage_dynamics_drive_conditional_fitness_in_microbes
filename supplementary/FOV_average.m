clear;
% name='pos1_2';
%
% postprocfolder = ['E:\CNphaC_005_001\output1_',name];
% tifseq = ['E:\CNphaC_005_001\',name];
% isTrechSet=0;

% name='pos2_1';
% postprocfolder = ['E:\',folder,'\output1_',name];

% folder='sample_from_f17_0dot5';
% folder='\both_dec_1to40_inc_41to80_4-2-1-0dot5-0dot1\dec3';
% thres_bright=1000;
% thres_dark=70;
% areathres = 300;

% folder='change_C_2_1_0dot5_0dot01_0dot1_0dot5_N_0dot5'

% folder='G:\summary_sim_change_Acetate_dec\summary_sim';

% folder = 'G:\both_inc_dec_2_to_0dot1_new\summary_sim_inc';
% incsign = 1;
% Lane = 0;

% folder = 'G:\control_change_same_media_20fovs\summary_sim_dec';
% Lane = 2;

% folder = 'G:\control_change_same_media_20fovs\summary_sim_inc';

% folder = 'D:\both_inc_dec_2_to_0dot1_new\dec\summary_sim'
Lane = 1;
% folder='D:\dec_fold_change_detection_lane1_3_2_1_0dot67_lane2_3_1dot5_1_0dot5_third\summary_sim';
% folder='F:\Experiment-346_change_ion\lane1\summary_sim';
% folder='E:\Experiment-299_Lane1_01_001_0001_Lane2_1_01_001-V7\lane1\summary_sim';
% folder='D:\control_change_same_media_20fovs\summary_sim_dec';

folder = 'F:\both_inc_dec_2_to_0dot1_new_cropped\inc\summary_sim';
% incsign = 1;
% Lane = 0;

% tic;
tic;
% folder='CNphaC_005_001';
postprocfolder=[folder,filesep];

isTrechSet = 1;


% position=1;
% if ~isTrechSet
%     load(fullfile(postprocfolder,sprintf('Position%06d.mat',position)));
% else
%     res=sim_readTrenches(postprocfolder);
% end
toc;


% start_time=1;
% end_time=900;
% smedia=[0,268,550,993];

%change 4 con
start_time=1;
end_time=830;
smedia=[0,253,463,692,961];
thres_bright = 300;
thres_dark = 100;
areathres = 400;
areathresup = 800;
minsize=0;

% thres_bright = 134; % 250, 200
% thres_dark = 34; % 100, 50



%control same media
% thres_bright = 1200;
% thres_dark = 300;
% areathres = 400;
% areathresup = 800;
% start_time=1;
% end_time=820;
% minsize = -1;
% Lane = 1;Lane = 2;


%change acetate from 2 to 1 to 0.5
% start_time=1;
% end_time=678;
% smedia=[0,140,394,678];
% Lane = -1;
% minsize = -1;
%
% thres_bright = 400;
% thres_dark = 80;
% areathres = 400;
% areathresup = 800;



%for sample 4
% thres_bright=280;
% thres_dark=100;
% areathres = 400;
% areathresup = 800;
% start_time=30;
% end_time=159;

% start_time=160;
% end_time=431;

% start_time=1;
% end_time=1231;
% smedia=[0,225,441,722,985,1231];






% start_time=1;
% end_time=372;
% smedia=[0,295,586,895,1161,1285];
tic;
drawind = 1;
pp=1;
thres_gen = 50;


sample_granules=cell(end_time-start_time+1,1);
sample_granules_pixel=cell(end_time-start_time+1,1);
sample_area=cell(end_time-start_time+1,1);
sample_rate_smooth=cell(end_time-start_time+1,1);
sample_rate=cell(end_time-start_time+1,1);
sample_divs=zeros(end_time-start_time+1,1);
sample_outside_gran_without=cell(end_time-start_time+1,1);
sample_outside_gran=cell(end_time-start_time+1,1);
sample_totf_outside_gran=cell(end_time-start_time+1,1);
sample_dilute=cell(end_time-start_time+1,1);
max_gran_size=10000;
rate_granule_size=cell(1,max_gran_size);
smooth_rate_granule_size=cell(1,max_gran_size);
darknum=0;
darknum1=0;
outn=0;
res1=[];


dd=dir(postprocfolder);
output=[];
chamber = 0;
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
            
            %             res1{chamber}.lineage{val}.daughters=lineage{val}.daughters;
            %             res1{chamber}.lineage{val}.length=lineage{val}.length;
            %             res1{chamber}.lineage{val}.framenbs=lineage{val}.framenbs;
            %
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
                ouside_gran=zeros(et-st+1,1);
                mgran=zeros(et-st+1,1);
                mean_outside_gran=zeros(et-st+1,1);
                totf_outside_gran=zeros(et-st+1,1);
                mean_outside_gran_dilute = zeros(et-st+1,1);
                for i=st:et
                    ia=ia+1;
                    if iscell(fluo)
                        gran(ia) = sum(fluo{i}>thres_bright);
                        mgran(ia) = max(fluo{i});
                        mean_outside_gran(ia)= mean(fluo{i}(fluo{i}<=thres_bright));
                        totf_outside_gran(ia)= sum(fluo{i}(fluo{i}<=thres_bright));
                    else
                        gran(ia) = sum(fluo(i,:)>thres_bright);
                        mgran(ia) =  max(fluo(i,:));
                        mean_outside_gran(ia) = mean(fluo(i,fluo(i,:)<=thres_bright));
                        totf_outside_gran(ia) = sum(fluo(i,fluo(i,:)<=thres_bright));
                    end
                end
                
                dau = lineage{val}.daughters(st:et);
                
                divs=find(dau);
                sf=1;
                
                for i=1:length(divs)
                    ef=divs(i);
                    
                    ff=localFramenbs(ef)-start_time+1;
                    sample_divs(ff) = sample_divs(ff)+1;
                end
                
                for i=1:length(gran)
                    if mgran(i)>thres_dark && area(i) > areathres
                        ff=localFramenbs(i)-start_time+1;
                        
                        out=[sample_granules{ff},gran(i)/area(i)];
                        sample_granules{ff}=out;
                        
                        out=[sample_area{ff},area(i)];
                        sample_area{ff}=out;
                        
                        out=[sample_granules_pixel{ff},gran(i)];
                        sample_granules_pixel{ff}=out;
                        
                        if mgran(i)<=thres_bright
                            out=[sample_outside_gran_without{ff},mean_outside_gran(i)];
                            sample_outside_gran_without{ff}=out;
                            
                            out=[sample_totf_outside_gran{ff},totf_outside_gran(i)];
                            sample_totf_outside_gran{ff}=out;
                        end
                        
                        out=[sample_outside_gran{ff},mean_outside_gran(i)];
                        sample_outside_gran{ff}=out;
                        
                    else
                        darknum=darknum+1;
                    end
                end
                
                outgran=[];
                for i=1:length(divs)
                    ef=divs(i);
                    dt=ef-sf+1;
                    if dt>thres_gen
                        break;
                    end
                    
                    if sf<ef-1
                        traj=gran(sf:ef-1);
                        straj=smooth((sf:ef-1)',traj,0.7);
                        outgran=[outgran;straj];
                        
                        
                        for ti=sf:ef-1
                            ff=localFramenbs(ti)-start_time+1;
                            
                            tti=ti-sf+1;
                            if mgran(ti)>thres_dark && area(ti) > areathres
                                
                                if tti>1 && dau(ti)==0 && straj(tti)>minsize %&& log(area(i)/area(i-1))*12>0.03
                                    out=[sample_rate_smooth{ff},straj(tti)-straj(tti-1)];
                                    sample_rate_smooth{ff}=out;
                                    
                                    out=[sample_rate{ff},traj(tti)-traj(tti-1)];
                                    sample_rate{ff}=out;
                                end
                                
                                if tti>1 && dau(ti)==0 %&& straj(tti)>=straj(tti-1)  %&& log(area(i)/area(i-1))*12>0.03
                                    
                                    gsize=int32(traj(tti-1)+1);
                                    out=[rate_granule_size{gsize},traj(tti)-traj(tti-1)];
                                    rate_granule_size{gsize}=out;
                                    
                                    gsize=int32(straj(tti-1)+1);
                                    out=[smooth_rate_granule_size{gsize},straj(tti)-straj(tti-1)];
                                    smooth_rate_granule_size{gsize}=out;
                                    %                                 gsize=int32(straj(tti)+1);
                                    %                                 out=[rate_granule_size{gsize},straj(tti)-straj(tti-1)];
                                    
                                end
                                
                                %                                 if ti>1 && mgran(i)<=thres_bright
                                %
                                %
                                %                                     out=[sample_dilute{ff},log(area(ti)/area(ti-1))*mean_outside_gran(ti)];
                                %                                     sample_dilute{ff}=out;
                                %                                 end
                                
                                
                            else
                                darknum1=darknum1+1;
                            end
                        end
                        
                    else
                        outgran=[outgran;nan];
                    end
                    sf=ef;
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
                    
                    plot(localFramenbs,mean_outside_gran,'-','Color',[0.9290 0.6940 0.1250]);
                    
                    set(gca,'FontName','Times new roman','FontSize',15)
                    title(['pos=',num2str(1),' trench=',num2str(chamber),' val=',num2str(val),' # div=',num2str(sum(dau>0))])
                    box on;
                end
            end
            
        end
    end
end
display(darknum);
display(darknum1);
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

statpromo=zeros(length(sample_outside_gran_without),2);
for i=1:length(sample_outside_gran_without)
    statpromo(i,1)=nanmean(sample_outside_gran_without{i});
    statpromo(i,2)=nanvar(sample_outside_gran_without{i});
end

statgran_pixel=zeros(length(sample_granules_pixel),2);
for i=1:length(sample_granules_pixel)
    statgran_pixel(i,1)=nanmean(sample_granules_pixel{i});
    statgran_pixel(i,2)=nanvar(sample_granules_pixel{i});
end

statpromo_totf=zeros(length(sample_totf_outside_gran),2);
for i=1:length(sample_totf_outside_gran)
    statpromo_totf(i,1)=nanmean(sample_totf_outside_gran{i});
    statpromo_totf(i,2)=nanvar(sample_totf_outside_gran{i});
end

%%

midstep=1;


meanval=zeros(max_gran_size/midstep,1);

range=2:100;
xx=range*midstep-1;

for i=1:(max_gran_size/midstep)
    
    tmp=[];
    for j=1:midstep
        indexnum=(i-1)*midstep+j;
        tmp=[tmp,rate_granule_size{indexnum}];
    end
    
    %     if i<2
    %         figure;
    %         [y,x]=hist(tmp,-20:0.1:20);
    %         bar(x,y/sum(y))
    %     end
    
    %     if i<30
    % %         min(tmp)
    %     end
    %     tmp=tmp(abs(tmp)<100);
    %     tmp=tmp( tmp>-5 & tmp<5);
    meanval(i)=mean(tmp);
end

% plot(xx,meanval(range),'o');

% plot(f,xx,meanval(range),'o')
% ylim([0 10])

%
% f = @(p,x) p(1).*x.^2+p(2).*x+p(3);
f = @(p,x) p(1).*(x+p(2)).^(-0.5)  + p(3);
options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-20,'MaxIter',1e20,'TolX',1e-20);
P = fminsearch(@(p) norm(meanval(range) - f(p,xx')),rand(3,1)*5,options)
P = [103.8371,16.0351,-11.123];
figure
% plot(xx',meanval(range), 'p')
hold on
plot(xx', f(P,xx'), '-k')


% plot(xx',300*(xx+10).^(-1.1))
% figure;
% f=fit(log(xx+5)',log(meanval(range)),'poly1');
% plot(f,log(xx)',log(meanval(range)),'o');

% plot(log(x(1:100))',log(meanval(1:100)),'o')
% hold on;
% plot(f,log(xx)',log(meanval(range)),'o')

% figure;
% f=fit(xx',meanval(range),'poly2')
% plot(f,xx',meanval(range),'o')
% hold on;
% xlabel('granule size')
% ylabel('mean production rate (rate>=0)')
% title('log-log plot')
% figure;
% plot(xx,104*(xx+16).^(-1/2)-11.11,'k');
% xlabel('granule size')
% ylabel('mean production rate (rate>=0)')
%
% figure;
% plot(xx,smooth(meanval(xx),0.4),'-',xx,0.2725*(xx).^(2/3))

% plot(xx,meanval(range),'o',xx,0.55*(xx).^(1/2))
%
% hold on;
% plot(xx,meanval(range),'o',xx,0.27*(xx).^(2/3))

% %%
% frames_plot=[200,260,300];
% 
% outplot=[];
% for i=1:length(frames_plot)
%     
%     val=frames_plot(i);
%     
%     
%     %     histv=sample_granules{val};
%     %     [y,x]=hist(histv(histv>0),0:0.02:0.8);
%     
%     histv=sample_rate_smooth{val};
%     [y,x]=hist(histv(abs(histv)<=20),-20:2:20);
%     
%     outplot=[outplot;y/sum(y)];
%     % plot(sample_granules{val}./sample_area{val},sample_area{val},'x','Markersize',2);
%     % ylim([0 6000])
%     % xlim([0 1.1])
%     aa=histv(histv>0);
%     sqrt(var(aa))/mean(aa)
% end
% 
% figure;
% leg=cell(size(frames_plot));
% 
% for i=1:length(frames_plot)
%     leg{i}=sprintf('Frame:%d',frames_plot(i));
%     
% end
% 
% 
% 
% 
% bar(x,outplot)
% legend(leg)
% 
% ylabel('Probability')
% xlabel('PHB granule over cell area')
% % xlabel('Production rate (# pixel per 5min)')
% 
% set(gca,'FontName','Times new roman','FontSize',20)
% 
% 
% 
% ylim([0 0.4])

%%

frames_plot=[300,400,500];
outplot=[];

k=-1;
for i=1:length(frames_plot)
    
    val=frames_plot(i);
    
    
    histv=sample_rate_smooth{val};
    [y,x]=hist(histv(abs(histv)>k),-20:1:20);
    
    histv=sample_rate{val};
    [y0,x0]=hist(histv(abs(histv)>k),-20:1:20);
    
    outplot=[outplot;y/sum(y)];
    % plot(sample_granules{val}./sample_area{val},sample_area{val},'x','Markersize',2);
    % ylim([0 6000])
    % xlim([0 1.1])
    
end


figure;
leg=cell(size(frames_plot));

for i=1:length(frames_plot)
    leg{i}=sprintf('Frame:%d',frames_plot(i));
    
end




bar(x,outplot)
legend(leg)
ylim([0 0.4])

ylabel('Probability')
xlabel('Production rate (# pixel per 10min)')

set(gca,'FontName','Times new roman','FontSize',20)
% title('Renormalized after removing 0 rates')

mean_num = sum(numcell) / 830;

%%
figure;

yyaxis left
% plot(start_time:end_time,stat(:,1),'-','linewidth',1)
hold on
plot(start_time:end_time,statpromo(:,1),'-','linewidth',1)
% title('Decreasing from 2g/L to 0.5g/L Acetate');
% title('Decreasing from 2g/L to 0.1g/L NH4Cl');
% title('Lane 2');


% if Lane==1
%     title('Lane 1');
%     con=[1,1,0.5,1];%lane 1
%     smedia=[0,254,396,686,820];%lane 1
% else
%     if Lane==2
%         title('Lane 2');
%         con=[1,0.5,1,1]; %lane 2
%         smedia=[0,254,541,686,820]; %lane 2
%
%     else
%         if Lane==0
%             if incsign==0
%                 con=[2,1,0.5,0.1];
%             else
%                 con=[0.1,0.5,1,2];
%             end
%         else
%             con=[2,1,0.5];
%         end
%
%     end
% end

% title('Increasing from 0.1g/L to 2g/L NH4Cl');
hold on
% plot(start_time:end_time,sqrt(stat(:,2)),'k-','linewidth',1)

set(gca,'FontName','Times new roman','FontSize',15,'linewidth',1)
ylim([0,0.20])
% ylim([0,0.3])
% yticks(0:0.2:0.6);
% yticklabels({'0%','20%','40%','60%'});


xlim([start_time,end_time]);
xlabel('Time (frame/10min)')

yyaxis right

plot(start_time:end_time,numcell,'linewidth',1)

plot(start_time:end_time,smooth(sample_divs,0.1),'r-','LineWidth',2)
plot(start_time:end_time,sample_divs,'r*')
%% 

figure;
divpercells=sample_divs/numcell;
plot(start_time:end_time,divpercells,'linewidth',1)

% ylim([0,5000])
set(gca,'FontName','Times new roman','FontSize',15)


% con=[4,2,1,0.5,0.1];
%
%
% con=[2,1,0.5,0.1,0.01];

% con=[2,1,0.5,0.1,0.01];

con=[2,1,0.5,0.1,0.01];
%
for i=1:length(smedia)-1
    xline(smedia(i),'linewidth',1)
    txt = sprintf("%.1fg/L",con(i));
    text(smedia(i),1200,txt,'HorizontalAlignment','left')
end

legend('Mean','Standard deviation','# of cells')


% figure;
% plot(start_time:end_time,sqrt(stat(:,2))./stat(:,1),'linewidth',2)
% title('CV');
%
% ylim([0,4])
%
% xlabel('Time (frame/5min)')
% xlim([start_time,end_time]);
% set(gca,'FontName','Times new roman','FontSize',15)
%
% %title([sprintf('(%d to %d)',start_time,end_time),'  mean=',sprintf('%.2f',nanmean(dlen)*0.055),'  CV=',sprintf('%.2f',sqrt(nanvar(dlen))/nanmean(dlen)),'  #sample=',num2str(sum(~isnan(dt)))])
%
%
% %%
% figure;

% plot(start_time:end_time,statpromo(:,1),'linewidth',2)
% ylim([120 220])
set(gca,'FontName','Times new roman','FontSize',15)

% plot(start_time:end_time,statgran_pixel(:,1),'m-','linewidth',2)
% plot(start_time:end_time,statarea(:,1)*0.001.*statpromo(:,1),'k-','linewidth',2)
% plot(start_time:end_time,statarea(:,1),'m-','linewidth',2)
% title('Cell area');
hold on;
% plot(start_time:end_time,sqrt(statarea(:,2)),'k-','linewidth',2)
% legend('Mean','Standard deviation')

% ylim([0 200])
ylabel('# of pixels')
xlim([start_time,end_time]);
xlabel('Time (10min per frame)')

set(gca,'FontName','Times new roman','FontSize',15)
%
%
%
% figure;
% plot(start_time:end_time,sqrt(statarea(:,2))./statarea(:,1),'linewidth',2)
% title('CV');
% set(gca,'FontName','Times new roman','FontSize',15,'linewidth',1)
% xlim([start_time,end_time]);
% ylim([0 1])
% xlabel('Time (frame/5min)')
%
%
%% 
figure;
plot(start_time:end_time,statpromo(:,1),'-','linewidth',1)
%%
% clear res;
save dec_less_smooth.mat;