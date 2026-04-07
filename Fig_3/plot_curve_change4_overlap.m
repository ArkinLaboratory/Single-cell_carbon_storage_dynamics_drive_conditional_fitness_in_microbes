clear;
tic;

%
files={'change4_dec_20FOVs.mat','change4_inc_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};
leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,0,0]);
tit={'Decreasing NH_4Cl','Increasing NH_4Cl'}
files=files(sel)
leg=leg(sel);

k=-1;

Lane=0;
figure;
sel=[3,2,1];

leg{1}={'2g/L to 1g/L','1g/L to 0.5g/L','0.5g/L to 0.1g/L'};
leg{2}={'0.1g/L to 0.5g/L','0.5g/L to 1g/L','1g/L to 2g/L'};

for f=1:length(files)
    load(files{f},'sample_divs','numcell','statarea','stat','start_time','end_time','con','smedia');
    %change 4 con
    % start_time=1;
    % end_time=961;
    % smedia=[0,253,463,692,961];
    % thres_bright = 300;
    % thres_dark = 80;
    % areathres = 400;
    % areathresup = 800;
    % minsize=-1;
    
    
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
    % leg={'2g/L to 1g/L','1g/L to 0.5g/L','0.5g/L to 0.1g/L'};
    
    %  leg={'2g/L to 1g/L','0.1g/L to 0.5g/L'};
    %  leg={'2g/L to 1g/L','1g/L to 2g/L'};

    %%
    
    
    
    % dd=smedia;
    % dd=dd(2:end-1);
    % dd=[23,dd];
    % numcell(dd+1)=numcell(dd);
    % numcell(dd+2)=(numcell(dd)+numcell(dd+4))/2;
    % numcell(dd+3)=numcell(dd+4);
    
    % stat(dd+1,:)=stat(dd,:);
    % stat(dd+2,:)=(stat(dd,:)+stat(dd+4,:))/2;
    % stat(dd+3,:)=stat(dd+4,:);
    
    % yyaxis left
    figure;
    hold on;
    head=0;
    for i=1:3
        range=(smedia(i+1)+head):(smedia(i+1)+200-1);
        plot((1:200-head)/6,stat(range,1),'linewidth',2);
        % title('Decreasing from 2g/L to 0.5g/L Acetate');
        % title('Decreasing from 2g/L to 0.1g/L NH4Cl');
        % title('Lane 2');
    end
    plot((1:200)*10/60, 0.059*ones(200-head,1),'k--','linewidth',1)
    
    
    xlim([0,(200-head)/6]);
    ylim([0.04,0.08]);
    yticks([0.04,0.06,0.08]);
    yticklabels({'4%','6%','8%'});
    xlabel('Time (hours)')
    set(gca,'FontName','Times new roman','FontSize',30,'linewidth',1)
    box off;
    % title('Increasing from 0.1g/L to 2g/L NH_4Cl');
    % hold on
    % plot(start_time:end_time,sqrt(stat(:,2)),'k-','linewidth',1)
    %
    % set(gca,'FontName','Times new roman','FontSize',15,'linewidth',2)
    % ylim([0.0,0.1])
    % ylim([0,0.3])
    % yticks(0:0.2:0.6);
    % yticklabels({'0%','20%','40%','60%'});
    
    
    
    % yyaxis right
    
    
    % plot(start_time:end_time,numcell,'linewidth',1)
    
    % plot(start_time:end_time,smooth(sample_divs,0.1),'r-','LineWidth',2)
    % plot(start_time:end_time,sample_divs,'r*')
    % ylim([0,10000])
    % set(gca,'FontName','Times new roman','FontSize',20)
    
    
    
    title(tit{f})
    % title('Comparison')
    
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
    % plot(start_time:end_time,statarea(:,1),'m-','linewidth',2)
    
    % title('Cell area');
    % hold on;
    % plot(start_time:end_time,sqrt(statarea(:,2)),'k-','linewidth',2)
    % legend('Mean','Standard deviation')
    
    % ylim([0 1500])
    % plot(start_time:end_time,smooth(sample_divs,0.1),'r-','LineWidth',2)
    % plot(start_time:end_time,sample_divs,'r*')
    % for i=1:length(smedia)-1
    %     xline(smedia(i),'linewidth',1)
    %     txt = sprintf("%.1fg/L",con(i));
    %     text(smedia(i),600,txt,'HorizontalAlignment','left','FontName','Times new roman','FontSize',20)
    % end
    
    
    legend(leg{f})
end
