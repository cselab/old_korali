% function [y, theta] = make_data(N, x, FLAG)

folder_name = './';

addpath('../../models/logistic/');

FLAG = 1;

N = 10;
x = repmat( 0:0.5:10, N,1 );


Nx = size(x,2);

y = zeros(N,Nx);
theta = zeros(N,3);

for i = 1:N
   
%     r = binornd(1,0.5);
%     theta(i,1) = r*300 + (1-r)*150;
    
    theta(i,1) = normrnd(200,20);
    theta(i,2) = normrnd(40,10);
    theta(i,3) = normrnd(1,0.1);
%     theta(i,3) = lognrnd(1,0.2);
    theta(i,4) = 5;
    
    y(i,:) = my_model( x(i,:), theta(i,1:3), FLAG );
    
    
    error = normrnd( 0, theta(i,4), 1,Nx);
    y(i,:) = y(i,:) + error;

    data.x = x(i,:);
    data.y = y(i,:);
    data.theta = theta(i,1:3);
    data.std_data = theta(i,4);
    
    file_name = [folder_name 'data_set_' sprintf('%03d',i) '.mat'];
    save(file_name, 'data');
    
    
end


all_data = [ reshape( repmat(1:N,Nx,1),N*Nx,1 ) , reshape( x',N*Nx,1 ) , reshape( y',N*Nx,1 )];

fileID = fopen('all_data.txt','w');
fprintf(fileID,'%6s \t %6s \t %s\n','ID','time','y');
fprintf(fileID,'%f \t %f \t %f\n', all_data');
fclose(fileID);

clf
plot(x',y','-o','LineWidth',2)
ax=gca;

ax.FontSize = 15;


grid on