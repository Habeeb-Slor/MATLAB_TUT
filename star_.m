n =5;

for i = 1:n
    line = '';
    for j = 1:i
        line = [line '*'];
    end
    disp(line)
end