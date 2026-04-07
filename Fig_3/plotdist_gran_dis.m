

clear;
tic;

files={'change4_dec_20FOVs.mat','change4_dec_20FOVs.mat','change4_dec_new.mat','change4_inc_20FOVs.mat','change4_inc_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};
% leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,0,0,0,0]);
tit={'Decreasing from 2g/L to 0.1g/L NH_4Cl','Increasing from 0.1g/L to 2g/L NH_4Cl'};
files=files(sel);


outplotnew=[];
outsample={};


thres=-1;
% smedia=[0,253,463,692,961];
for f=1:length(files)
    load(files{f},'sample_granules','stat');
    % framesplot=[250,461];
    framesplot=[250,320];
    % framesplot=[250,320];
    
%     for i=1:length(framesplot)
        val=framesplot(f);
        histv=sample_granules{val};
        % histv=histv(histv>0);
%         sqrt(var(histv))/mean(histv)
        
        mean(histv);
        [y,x]=hist(histv,0:0.05:0.7);
        outplotnew=[outplotnew;y/sum(y)];
        outsample{f}=histv;
%         
%         outsample{i}=histv(abs(histv)>thres & abs(histv)<=20);
       
%         outputmean(i)=mean(histv(abs(histv)>thres & abs(histv)<=20));
%     end
end



leg=cell(size(framesplot));

% ss=1:2;
% for i=1:length(framesplot)
%     leg{i}=sprintf('Frame:%d (before %d switch)',framesplot(i),ss(i));
%     
% end


figure('Renderer', 'painters', 'Position', [10 20 450 250]);
bar(x,outplotnew)
% legend(leg)
% ylim([0 0.5])
% xlim([0 1])
% yyaxis left
% plot(1:961,stat(:,1),'b');
% 
% yyaxis right
% plot(1:961,outputmean);

% ylabel('Probability')
% xlabel('PHB granule size/cell size')
% yticks([])
% xticks([])
% title(tit{1})
% xlim([0,0.8])
set(gca,'FontName','Times New Roman','FontSize',20,'Linewidth',1)
ylabel('Probability')
xlabel('PHB granule size over cell size')
box off
% set(gca,'FontName','Times new roman','FontSize',20,'linewidth',2)
toc;


[h,p]=kstest2(outsample{1},outsample{2});

[out1,b]=hist(outsample{1},0:0.01:1);
[out2,b]=hist(outsample{2},0:0.01:1);
out1=out1./sum(out1);
out2=out2./sum(out2);
ret=bhattacharyya(out1,out2)
p
ret1=bhattacharyyaDistance(out1,out2)
% mean(outsample{1})
% mean(outsample{2})
% mean(outsample{3})
% [h,p]=kstest2(outsample{1},outsample{2});p
% [h,p]=kstest2(outsample{1},outsample{3});p