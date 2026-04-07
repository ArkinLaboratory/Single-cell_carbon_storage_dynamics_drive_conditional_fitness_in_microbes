clear;
tic;

%
% files={'change4_dec_20FOVs.mat','change4_inc_20FOVs.mat','control_lane2_new.mat','control_lane1_new.mat'};

files={'change4_dec_new.mat','change4_inc_new.mat','control_lane2_new.mat','control_lane1_new.mat'};
leg={'(Exp.1) 2g/L NH_4Cl','(Exp.1) 0.1g/L NH_4Cl','(Exp.2 lane 2) 1g/L NH_4Cl','(Exp.2 lane 1) 1g/L NH_4Cl'};
sel=logical([1,1,1,1]);

files=files(sel)
leg=leg(sel);

k=-1;
start_time=180;
end_time=240;
dt=[];
dlen=[];
outplotnew=[];
outplotnew1=[];
minframes=20;
pp=1;


for f=1:length(files)
    load(files{f},'res1');
    outnum=0;
    res=res1;
    for chamber=1:length(res)
        if ~isempty(res{chamber})
            for val=1:length(res{chamber}.lineage)
                allframe=length(res{chamber}.lineage{val}.framenbs);
                framenbs = res{chamber}.lineage{val}.framenbs;
                st=find(framenbs>=start_time,1);
                et=find(framenbs>end_time,1);
                
                if isempty(et)
                    et=allframe;
                else
                    et=et-1;
                end
                
                
                if ~isempty(st) && et-st+1>=minframes
                    
                    y=res{chamber}.lineage{val}.length(st:et);
                    dau=res{chamber}.lineage{val}.daughters(st:et);
                    
                    
                    if pp
                        for t=st+1:length(y)
                            if (y(t)-y(t-1))/y(t-1)>0.4 && t<length(y) && dau(t+1)>0
                                y(t)=y(t-1);
                                %y(t)=y(t-1)+abs(y(t-1)-y(t-2));
                                if t<length(y)
                                    dau(t+1)=0;
                                end
                            end
                            
                        end
                        
                        
                    end
                    
                    if sum(dau>0)>1
                        
                        times=find(dau>0);
                        dt=[dt,diff(times)];
                        
                        if times(1)==1
                            dlen=[dlen,y(times(2:end)-1)];
                        else
                            dlen=[dlen,y(times-1)];
                        end
                        outnum=outnum+1;
                        
                    end
                    
                end
            end
        end
    end
    dt(dt<=4)=nan;
    dt(dt>50)=nan;
    
    dlen(dlen*0.055>10)=nan;
    
    [y,x]=hist(dt,0:2:60);
    [y1,x1]=hist(dlen*0.055,0:0.3:10);
    
    
    display([sprintf('(%d to %d)',start_time,end_time),'  mean=',sprintf('%.2f',nanmean(dt)),'  CV=',sprintf('%.2f',sqrt(nanvar(dt))/nanmean(dt)),'  #sample=',num2str(sum(~isnan(dt)))])
    display([sprintf('(%d to %d)',start_time,end_time),'  mean=',sprintf('%.2f',nanmean(dlen)*0.055),'  CV=',sprintf('%.2f',sqrt(nanvar(dlen))/nanmean(dlen)),'  #sample=',num2str(sum(~isnan(dlen)))])

    outsample{f}=dt;
    outsample1{f}=dlen;
    outplotnew=[outplotnew;y/sum(y)];
    outplotnew1=[outplotnew1;y1/sum(y1)];
end


%%
figure('Renderer', 'painters', 'Position', [10 20 500 200]);
  
bar(x/6,outplotnew)
% legend(leg)
ylim([0 0.4])
xlim([1-0.2,4.2]);
box off;
ylabel('Probability')
xlabel('Doubling time (hours)')
% title(['Frame: ',num2str(start_time),' to ',num2str(end_time)])
set(gca,'FontName','Times New Roman','FontSize',26.8,'Linewidth',1)
toc;

figure('Renderer', 'painters', 'Position', [10 20 500 200]);

bar(x1,outplotnew1)
% legend(leg)
xlabel('Cell length at division (\mum)')
xlim([2-0.2,5.2]);
ylim([0 0.4])
box off;
ylabel('Probability')
% title(['Frame: ',num2str(start_time),' to ',num2str(end_time)])
set(gca,'FontName','Times New Roman','FontSize',26.8,'Linewidth',1)

[h,p]=kstest2(outsample1{1},outsample1{2})