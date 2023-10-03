tic
max = 2000000;
list = 1:max;

num = 2;
ind = 2;
while num < max
    list( mod(list,num) == 0 & list ~= num) = [];
    
    ind = ind + 1;
    if ind <= length(list)
        num = list(ind);
    else
        num = max;
    end
    
end

sum(list)
toc