clear

for i=1:5
    
    file_name = sprintf('../data/theta/theta_%03d.txt',i);
    d = load( file_name );
    
    plotmatrix_hist( d(:,1:end-2) );
   
    file_name = sprintf('theta_%03d.eps',i);
    saveas(gcf, file_name,'epsc')

end