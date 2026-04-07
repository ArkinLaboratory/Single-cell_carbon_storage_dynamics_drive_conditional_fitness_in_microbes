clear;
tic;

%
files={'change4_dec_20FOVs.mat','change4_inc_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};
leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,1,1]);
% tit={'Decreasing NH_4Cl','Increasing NH_4Cl'}
files=files(sel);
% leg=leg(sel);

% k=-1;

% Lane=0;
figure;
hold on;
% sel=[3,2,1];

% leg{1}={'2g/L to 1g/L','1g/L to 0.5g/L','0.5g/L to 0.1g/L'};
% leg{2}={'0.1g/L to 0.5g/L','0.5g/L to 1g/L','1g/L to 2g/L'};

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
%     head=0;
%     for i=1:3
%         range=(smedia(i+1)+head):(smedia(i+1)+200-1);
%         plot((1:200-head)/6,stat(range,1),'linewidth',2);
%         % title('Decreasing from 2g/L to 0.5g/L Acetate');
%         % title('Decreasing from 2g/L to 0.1g/L NH4Cl');
%         % title('Lane 2');
%     end
    

    confidence_level=0.95;
    alpha = 1 - confidence_level;
    
    range=1:250;
    
    for dd=1:length(numcell(range))
        t_score(dd) = tinv(1 - alpha/2, numcell(dd)-1);
    end
    
    
%     plot((range)*10/60,smooth(stat(range,1)+t_score'.*sqrt(stat(range,2))./sqrt(numcell(range)),0.01),'k-','linewidth',1)
%     plot((range)*10/60,smooth(stat(range,1)-t_score'.*sqrt(stat(range,2))./sqrt(numcell(range)),0.01),'k-','linewidth',1)
    plot((1:250)/6,stat(1:250,1),'linewidth',2);
    
    
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
    
    
    
 
  
    
    
    
end
plot((1:250)*10/60, 0.058*ones(250,1),'k--','linewidth',1)
xlim([0,250/6]);
ylim([0,0.20]);
yticks([0,0.05,0.1,0.15,0.2]);
yticklabels({'0%','5%','10%','15%','20%'});
xlabel('Time (hours)')
ylabel('Mean PHB granule area over cell area')
set(gca,'FontName','Times new roman','FontSize',20,'linewidth',1)
box off;
legend(leg)
