
close all;

for i=1:5
    
    c = [ 2 3 1 0.1] + 0.4*randn(1,4);
    c(4) = 0.3;
    
    disp(c)
    
    t = 0:0.1:2;
    y =  c(1)*sin( c(2)*t +c(3) ) + normrnd(0,c(4),1,length(t));

    p = plot(t,y,'-o');
    p.MarkerFaceColor = p.Color;
    p.MarkerSize = 6;
    grid on
    hold on
    
    ax = gca;
    ax.FontSize = 16;
    
    file_name = sprintf('%03d.dat',i);
    dlmwrite( file_name, [t' y'], '\t');
    
end

saveas(gcf,'data.eps','epsc')
