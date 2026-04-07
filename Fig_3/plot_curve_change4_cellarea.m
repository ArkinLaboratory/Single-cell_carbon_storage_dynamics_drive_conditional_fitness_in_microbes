clear;
tic;

%
files={'change4_inc_20FOVs.mat','change4_dec_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};
leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,0,0]);
tit={'Decreasing from 2g/L to 0.1g/L NH_4Cl','Increasing from 0.1g/L to 2g/L NH_4Cl'}
files=files(sel)
leg=leg(sel);

k=-1;

Lane=0;

pos1 = [0.1, 0.4, 0.4, 0.25];
pos2 = [0.1, 0.25, 0.4, 0.1];



for f=1:length(files)
    f
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
%     leg={'FOV mean of granule size/cell area','Standard deviation','# of cells'};
    leg={'Mean granule area/cell area','Smoothed curve','# divisions','Mean cell area'};

    fig=figure('Renderer', 'painters', 'Position', [10 20 1000 1000]);
    
    subplot1 = axes('Parent', fig, 'Position', pos1);
    yyaxis left
    hold on;
    
    confidence_level=0.95;
    alpha = 1 - confidence_level;
    
    for dd=1:length(numcell)
        t_score(dd) = tinv(1 - alpha/2, numcell(dd)-1);
    end

    col=[0.133, 0.545, 0.133]

    plot(subplot1,(start_time:end_time)*10/60,smooth(sample_divs./numcell,0.05),'b','Color', col,'linewidth',2)
    plot(subplot1,(start_time:end_time)*10/60,sample_divs./numcell,'.','Color', col)

    
    plot(subplot1,(start_time:end_time)*10/60,smooth(stat(:,1)+t_score'.*sqrt(stat(:,2))./sqrt(numcell),0.01),'-','linewidth',1)
    plot(subplot1,(start_time:end_time)*10/60,smooth(stat(:,1)-t_score'.*sqrt(stat(:,2))./sqrt(numcell),0.01),'-','linewidth',1)
   
    plot(subplot1,(start_time:end_time)*10/60,stat(:,1),'-','linewidth',1)


  

    ylim([0.01,0.15])
    yticks([])
    yticks(0:0.05:0.15);
    yticklabels({'0%','5%','10%','15%'});
    xlim([start_time,end_time]*10/60);
    xticks([])
    
    box off;

    yyaxis right
    hold on;
    
    % plot(subplot1,(start_time:end_time)*10/60,smooth(sample_divs,0.05),'r-',(start_time:end_time)*10/60,sample_divs,'r.','LineWidth',2)
    set(gca,'FontName','Times new roman','FontSize',20,'linewidth',1)

       
    plot(subplot1,(start_time:end_time)*10/60,smooth(statarea(:,1),0.01)+t_score'.*sqrt(statarea(:,2))./sqrt(numcell),'-','linewidth',1)
    plot(subplot1,(start_time:end_time)*10/60,smooth(statarea(:,1),0.01)-t_score'.*sqrt(statarea(:,2))./sqrt(numcell),'-','linewidth',1)
    
    plot(subplot1,(start_time:end_time)*10/60,statarea(:,1),'-','linewidth',1)
 
    
    ylim([250 1200])
    yticks(250:250:1000);
    
%     legend(leg)
    for i=1:length(smedia)-1
        xline(smedia(i)*10/60,'linewidth',1)
%         txt = sprintf("%.1fg/L",con(i));
%         text(smedia(i),0.04,txt,'HorizontalAlignment','left','FontName','Times new roman','FontSize',20)
    end
    
    
    subplot2 = axes('Parent', fig, 'Position', pos2);
    
    
    outcurve=ones(1,smedia(end));
    con0={[0.1,0.5,1.0,2.0],[2.0,1.0,0.5,0.1]};
    con=con0{f};
    for i=2:length(smedia)
        outcurve(smedia(i-1)+1:smedia(i))=outcurve(smedia(i-1)+1:smedia(i))*con(i-1);
    end
        
    plot(subplot2,(start_time:end_time)*10/60,outcurve,'linewidth',2)
    
    
    
    xlim([start_time,end_time]*10/60);
    xticks([0,40:40:160])
    box off;
    xlabel('Time (hours)')
    yticks([0.1,1.0,2.0]);
    yticklabels({'0.1','1.0','2.0'});
    ylabel('NH_4Cl (g/L)')
    set(gca,'FontName','Times new roman','FontSize',20,'linewidth',1)
end