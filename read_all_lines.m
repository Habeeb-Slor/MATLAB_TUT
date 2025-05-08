clear;
clc;

fileID = fopen("untitled10.txt");

i = 0;

while ~feof(fileID)

    

    if i == 1
        fgetl(fileID);
    else
        line = fgetl(fileID);
        disp(line);
    end
    
    i = i+1;
end