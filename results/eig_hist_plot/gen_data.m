nbins = 16;
[count, binvalues] = hist(data(:,1), nbins);

hold on
%bar(binvalues, count./sum(count), 'r', 'barwidth', 1);
bar(binvalues, count./sum(count), 'b', 'barwidth', 1);

count_ms = count./sum(count);
for i=1:length(count_ms)
    if (i ~= 3)
        count_ms(i) = 0;
    end
end
%bar(binvalues, count_ms, 'w', 'barwidth', 1);
bar(binvalues, count_ms, 'r', 'barwidth', 1);

title('Distribution of automatically designed categorization experiments');
ylabel('Frequency [%]\rightarrow');
xlabel('Expected information gain (E_{y|x}[D_{KL}]) [nat] \rightarrow');
set(gcf, 'Position', [100, 100, 400, 300]);
axis([0 0.3 0 0.25]);
box on;

%hold on
%bar(data(:,1));
%title('Information gain of automatically designed categorization experiments');
%ylabel('Expected information gain (E_y[D_{KL}]) [nat] \rightarrow');
%xlabel('Sorted designed experiment \rightarrow');
%set(gcf, 'Position', [100, 100, 400, 300]);
%axis([0 36115 0 0.35]);
%box on;

