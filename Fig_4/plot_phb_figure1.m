clear;
tic;

%
files={'change4_dec_20FOVs.mat','change4_inc_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};
leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,0,0]);
tit={'Decreasing from 2g/L to 0.1g/L NH_4Cl','Increasing from 0.1g/L to 2g/L NH_4Cl'}
files=files(sel)
leg=leg(sel);

k=-1;

Lane=0;

for f=1:length(files)
    f
    load(files{f},'sample_divs','numcell','statarea','stat','start_time','end_time','con','smedia');
    
    
%     leg={'Mean granule size over cell area'};
    
    figure;
    plot((start_time:end_time)*10/60,stat(:,1),'-','linewidth',2)
    plot((start_time:end_time)*10/60,stat(:,1)+sqrt(stat(:,2))./sqrt(numcell),'-','linewidth',2)
    
    
    hold on
    load(files{f},'sample_rate_smooth');
    set(gca,'FontName','Times new roman','FontSize',20)
    
    title(tit{f})
    
    for i=1:length(smedia)-1
        xline(smedia(i+1)*10/60,'linewidth',1)
        txt = sprintf("%.1fg/L",con(i));
        text(smedia(i+1)*10/60,0.1,txt,'HorizontalAlignment','right','FontName','Times new roman','FontSize',20)
    end
    
    
%     legend(leg)
    set(gca,'FontName','Times new roman','FontSize',20,'linewidth',2)
    ylim([0.00,0.15])
    
    yticks(0:0.05:0.15);
    yticklabels({'0%','5%','10%','15%'});
    xlim([start_time,end_time]*10/60);
    xlabel('Time (hours)')
    ylabel('Mean granule size over cell area')
    box off;
    
end