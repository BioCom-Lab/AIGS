function [] = fig_class(labels)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
n_class = max(labels);
for ii = 1:n_class
    a = find(labels==ii,1);
    b = a;
    c = length(find(labels==ii));
    d = c;
    rectangle('Position',[a-0.5,b-0.5,c,d],'EdgeColor',[0.5 0.5 0.5],'LineWidth', 2)
end
end

