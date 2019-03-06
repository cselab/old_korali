clc; clear;

% gp_file  = '../data/gp1.dat';
% out_file = '../data/out1.dat';
gp_file  = '../../build/examples/data/gp1.dat';
out_file = '../../build/examples/data/out1.dat';

addpath('./shplot')

k = 1; 

yt = load(out_file); 

[~,ind]=sort(yt(:,k));
indt = abs(yt(:,end-1))<1e-12;
yt(indt,end-1) = 0;

% subplot(2,1,1)

% subplot(2,1,2)
% semilogy( yt(ind,k), yt(ind,end-1) ); grid on


% load train data set
t = readtable( gp_file, 'HeaderLines', 12, 'ReadVariableNames', false,'Delimiter',' ');
train = table2array(t(:,1:end));

fig = figure(); clf

grid on; hold on



opt = {'Color', [0.9290    0.6940    0.1250] };
H = shplot( yt(ind,k), yt(ind,end-2), sqrt(yt(ind,end-1)) , opt{:});
H.line.LineWidth = 3;

% p3 = scatter( yt(ind,k), yt(ind,end-2) , 100, '.' );
% p3.MarkerFaceColor = [0.9290    0.6940    0.1250];
% p3.MarkerFaceAlpha = .5;




p1 = plot(train(:,k+1),train(:,1),'kx');
p1.MarkerFaceColor = [0.9290    0.6940    0.1250];
p1.LineWidth = 3;
p1.MarkerSize = 10;



p4 = plot(  yt(ind,k), yt(ind,end), 'b.-' );
p4.MarkerSize = 10;


lg=legend([p1,H.line,H.patch,p4],'data','GP mean','GP one std','exact');
lg.Location = 'best';

set(gca,'FontSize',20)