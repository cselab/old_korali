a = load(['../data/001.dat']);


for i=1:5
    b = load(['../data/00' num2str(i) '.dat']);
    plot(a(:,1),b(:,2),'-o');
    hold on
end