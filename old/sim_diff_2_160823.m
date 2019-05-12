% Simulate the restircted / hindered diffusion and the exchange between
% them.

clear

main
%% set figure
f = figure(100); cla
a1 = subplot(3,1,1); hold on; hold on; set(a1,'Ycolor','w','Xcolor','w');
a2 = subplot(3,1,2); ta2='Displacement along x'; title(ta2); set(a2,'fontsize',10);
a3 = subplot(3,1,3); hold on; hold on; title('Apparent Diffusion Coefficient'); set(a3,'fontsize',10); xlabel('time (ms)')

set(a1,'position',[0.03 0.43, 0.94, 0.55],'xtick',[],'ytick',[]);
set(a2,'position',[0.07 0.27, 0.86 0.13]);
set(a3,'position',[0.07 0.04, 0.86 0.17]);

%% parameters for the molecules
axes(a1);cla;
axes(a2);cla;
axes(a3);cla;

N = 2000;       % number of molecules.
STEP = 300;    % number of repetition.
A = 1;       % mean displacement for 1 step (gaussian probability).
trans = 0.00;   % probability of penetrating.
theta = rand(N,1)*2*pi; % displacement directions for the first step.

ratioD = 0.0010; % the ratio of the unit length in this script to mm.
ratioT = 0.00012; % the ratio of the unit time (1 STEP) to sec.
%% parameters for the cells
% center of the cells.
cen = packing.final_positions{1}'; cen=cen(:,[2 1]); cen=cen-repmat([mean(cen(:,1)) mean(cen(:,2))],[size(cen,1) 1]);

% Radius of the "cells"
R = [axons.d{1}.*axons.g_ratio{1} axons.d{1}]; % internal / external diameter in um

t = linspace(0,2*pi);
subplot(a1);
for k =1:size(R,1)
    plot(R(k,1)*cos(t)+cen(k,1), R(k,1)*sin(t)+cen(k,2),'b','linewidth',1);
    plot(R(k,2)*cos(t)+cen(k,1), R(k,2)*sin(t)+cen(k,2),'b','linewidth',1);
end
axis equal

tp = 1.2 * max(abs(cen(:)));
% tp = tp*0.4;
% 
% xlim([-tp,tp]); 
% ylim([-tp,tp]); 

%% set startpoints

% set2: rectangle
yym=min(cen(:,2)); yyM=max(cen(:,2));
xxm=min(cen(:,1)); xxM=max(cen(:,1));

squareratio = .6;
Ntmp=0;
st1=zeros(N,1); st2 = st1;
while Ntmp<N
    st1(st1==0) =  (xxm+xxM)/2+squareratio*(rand(N-Ntmp,1)-.5)*(xxM-xxm); % x-coordinates of the starting points.
    st2(st2==0) =  (yym+yyM)/2+squareratio*(rand(N-Ntmp,1)-.5)*(yyM-yym); % y-coordinates
    for p = 1:N
        dist = ((cen(:,1)-st1(p)).^2 + (cen(:,2)-st2(p)).^2).^0.5;
        if any(dist>R(:,1) & dist<=R(:,2))
            st1(p)=0; st2(p)=0;
        end
    end
    Ntmp = nnz(st1);
end
% % set1: in a circle
% temp = rand(N,1); temp2 = rand(N,1);
% st1 = temp2.*cos(temp*2*pi)*R(1);
% st2 = temp2.*sin(temp*2*pi)*R(1);

% % set3: 
% st1 = 4 * (rand(N,1)-0.5)*2; % x-coordinates of the starting points.
% st2 = 4 * (rand(N,1)-0.5)*2; % y-coordinates
% st3 = zeros(N,1);  
% for p = 1:N
%     dist = ((cen(:,1)-st1(p)).^2 + (cen(:,2)-st2(p)).^2).^0.5 - R;
%     st3(p) = find(dist == min(dist),1).* max(dist<=0);
% end
% st1 = st1(st3==0); st2 = st2(st3==0);
% N = length(st1);

%%
fr(STEP+1) = struct('cdata',[],'colormap',[]);

%% draw initial points
st3 = zeros(N,1);   % intra or extra cell (state).
for p = 1:N
    dist = ((cen(:,1)-st1(p)).^2 + (cen(:,2)-st2(p)).^2).^0.5 - R(:,1);
    st3(p) = find(dist == min(dist),1).* max(dist<=0);
end

STATE = [st1,st2, st3, zeros(N,1)]; % x-coordinate, y-coordinate, state, phase

% keep only extraaxonal spins
STATE(STATE(:,3)==0,:)=[];
N = size(STATE,1);

INITIAL = STATE(:,1:3); %save the start parmeters (phases are zero).

axes(a1); 
p1 = scatter(STATE(INITIAL(:,3)>0,1),STATE(INITIAL(:,3)>0,2),'.r');
p2 = scatter(STATE(INITIAL(:,3)==0,1),STATE(INITIAL(:,3)==0,2),'.g');

fr(1) = getframe(f);
%% Simulation
tic
meandisp = zeros(1,STEP);
diffcoef = zeros(1,STEP);
savediff = zeros(N,STEP);
for s = 1:STEP
    flight = A * randn(N,1);
    theta = rand(N,1)*2*pi; % random direction
    for k = 1:N
        [res, cent, stat] = endpoint2_140805(R,cen,STATE(k,1:2),theta(k),flight(k),trans,STATE(k,3));
        STATE(k,1:2) = res;
        STATE(k,3) = stat;
    end
    
    figure(f)
    subplot(a1)
    delete(p1); delete(p2); 
    xlim([-tp,tp]); 
    ylim([-tp,tp]);
    axis equal
    p1 = scatter(STATE(INITIAL(:,3)>0,1),STATE(INITIAL(:,3)>0,2),'.r');
    xlim([-tp,tp]); 
    ylim([-tp,tp]); 
    axis equal
    p2 = scatter(STATE(INITIAL(:,3)==0,1),STATE(INITIAL(:,3)==0,2),'.g');
    xlim([-tp,tp]); 
    ylim([-tp,tp]); 
    axis equal
    
    savediff(:,s) = (((STATE(:,1)-INITIAL(:,1)).^2 + (STATE(:,2)-INITIAL(:,2)).^2)*ratioD^2)./(s*ratioT);
    disps=sqrt((STATE(:,1)-INITIAL(:,1)).^2 + (STATE(:,2)-INITIAL(:,2)).^2);
    meandisp(s) =mean(disps);
    diffcoef(s) = (meandisp(s)*ratioD)^2/(s*ratioT);
    subplot(a2);
    hold off
    histogram(STATE(:,1)-INITIAL(:,1),linspace(0,20,50)-10,'Normalization','pdf')
    t=title(ta2);
    xlim([-10 10]); ylim([0 .5]);
%    plot((1:s)*ratioT*1000,(meandisp(1:s)*ratioD).^2,'b-','linewidth',1.5); xlim([2*ratioT*1000,STEP*ratioT*1000*1.1]);
    subplot(a3);
    plot((1:s)*ratioT*1000,diffcoef(1:s),'r-','linewidth',1.5); xlim([2*ratioT*1000,STEP*ratioT*1000*1.1]); ylim([0 0.003])

    % commnet out next the lines to make the simulation faster.
    drawnow
    fr(s+1) = getframe(f);
end
toc
%% finish

D = ((meandisp(s) * ratioD)^2) / (STEP*ratioT);
fprintf('Diffusion coefficient: %d mm^2/sec\n',D);
   
    
    
    
    
    
    
    


