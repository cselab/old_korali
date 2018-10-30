% function [y, theta] = make_data(N, x, FLAG)

folder_name = './';

addpath('../../models/normal/');


N = 200;
Nx = 1;

x = repmat( 0 , N,1 );

y = zeros(N,1);
theta = zeros(N,2);

for i = 1:N
   

    theta(i,1) = normrnd(5,0.1);
    theta(i,2) = 1;
    
    y(i,:) = my_model( 0, theta(i,1) );
    
    
    error = normrnd( 0, theta(i,2), 1, Nx );
    y(i,:) = y(i,:) + error;

    data.x = 0;
    data.y = y(i,:);
    data.theta = theta(i,1:1);
    data.std_data = theta(i,2);
    
    file_name = [folder_name 'data_set_' sprintf('%03d',i) '.mat'];
    save(file_name, 'data');
    
    
end


all_data = [ reshape( repmat(1:N,Nx,1),N*Nx,1 ) , reshape( x',N*Nx,1 ) , reshape( y',N*Nx,1 )];

fileID = fopen('all_data.txt','w');
fprintf(fileID,'%6s \t %6s \t %s\n','ID','time','y');
fprintf(fileID,'%f \t %f \t %f\n', all_data');
fclose(fileID);

clf
plot(x',y')


