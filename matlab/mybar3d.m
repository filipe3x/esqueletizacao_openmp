%T = csvread('resultstimes16-12-2018-184658.csv');

fid = fopen('resultstimes16-12-2018-184658.csv','rt');
T = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',';');
fclose(fid);

fid = fopen('resultstimes18-12-2018-175111.csv','rt');
C = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',';');
fclose(fid);

P = cell2mat(C)./cell2mat(T);
P = P .* 100;
b = bar3(P);

xlabel('Procs'), ylabel('Image Sizes'), zlabel('communications ratio %');
set(gca,'XTick',0:2:32);
t = {'330','512','1024','1536','2048','4096','6144','8192'};
t = reshape (t, [8,1]);
set(gca, 'YTickLabel',t);