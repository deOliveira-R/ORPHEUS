batchSize = 100;
sampleSize=10000;
x = [ ];
for i = 1:sampleSize
    x = [x mean(randi(6,batchSize))];
end
f = figure('visible','off');
hist(x,32);
xlabel('Value');
ylabel('Frequency');
xlim([1 6])
saveas(f, 'Fig_01.pdf');


