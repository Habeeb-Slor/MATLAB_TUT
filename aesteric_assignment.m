clear;
clc;
n = input("enter number of rows: ");

for i = 1:n
    for j = 1:(n-i)
       fprintf(" ")
    end
    for k = 1 : i
        fprintf(" *")
    end
    fprintf("\n")
end