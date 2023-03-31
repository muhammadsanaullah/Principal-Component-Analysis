
load SpikeDataforPCA.mat

session(1).wf=double(session(1).wf);
session(2).wf=double(session(2).wf);

wf1 = session(1).wf; %waveforms from the 1st day
wf2 = session(2).wf; %waveforms from the 2nd day
stamps1 = session(1).stamps; %time stamps from the 1st day
stamps2 = session(2).stamps; %time stamps from the 2nd day

[coeff1,score1,latent1] = pca(wf1); %Compute principal components
[coeff2,score2,latent2] = pca(wf2);
edges{1}=[-300:5:300]; %Bin for the X-axis
edges{2}=[-250:5:250]; %Bin for the Y-axis
% the first 2 columns of score become X-Y axes, PC1 PC2
h=hist3(score1(:,1:2),edges); %Compute a 2-D histogram
s=surface(h'); %Visualize the histogram as a surface (note the apostrophe)
set(s,'XData',edges{1}) %Label the X-axis
set(s,'YData',edges{2}) %Label the Y-axis

temp_wf = score1; x = temp_wf(:,1); y = temp_wf(:,2);
euc_dist = sqrt((5 - x).^2 + (55 - y).^2);
ind = find(euc_dist < 70); % indices for which euc distance < threshold
% original space is found by score*eigenvectors
origSpace = score1(ind,:)*coeff2';
h2 = hist3(score1(ind,1:2),edges); figure;
s2=surface(h2'); %Visualize the histogram as a surface (note the apostrophe)
set(s2,'XData',edges{1}) %Label the X-axis
set(s2,'YData',edges{2}) %Label the Y-axis
% plotting average of those neuronal waveforms at the indices required
template = mean(origSpace);
figure; plot(template); 

% RMS calculation for spike-sorting
[r,c] = size(wf2); rms = []; tot = 0;
for i = 1:1:r
    for j = 1:1:c
        value = (template(j) - wf2(i,j))^2;
        tot = tot + value;
    end
    rms_val = sqrt(tot/r);
    rms = [rms rms_val];
    tot = 0;
end

% temp_s = h; tempx = edges{1,1}; tempy = edges{1,2};
% euc_dist = 0; [rx,cx] = size(tempx); [ry,cy] = size(tempy);
% for i = 1:1:cx
%     for j = 1:1:cy
%         euc_dist = sqrt((5-tempx(i))^2+(65-tempy(j))^2);
%         if (euc_dist > 70)
%             temp_s(i,j) = 0;
%         end
%     end
% end
% 
% figure; ts = surface(temp_s');
% set(ts,'XData',edges{1}) %Label the X-axis
% set(ts,'YData',edges{2}) %Label the Y-axis
% %figure; plot(mean(h)); hold on; plot(mean(temp_s)); 

data_retained = [];
for i = 1:1:size(latent1)
    temp = latent1(i) / (sum(latent1));
    data_retained = [data_retained temp];
end
data_retained = 1 - data_retained;
%data_retained = smooth(data_retained);
p = polyfit(1:1:48,data_retained,20);
y1 = polyval(p,1:1:48);
figure; plot(1:1:48,data_retained,'x');
%ylim([0 1]); xlim([0 48]);
hold on; plot(1:1:48,y1);
