clear;
clc;
n = 5;

for i = 1:5
    for j = 1:(5-i)
       fprintf(" ")
    end
    for k = 1 : i
        fprintf("*")
    end
    fprintf("\n")
end