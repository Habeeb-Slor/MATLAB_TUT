fileID = fopen('untitled10.txt', 'r');

if fileID == -1
    error('Cannot open file.');
end

while ~feof(fileID)
    ch = fread(fileID, 1, '*char');  % Read one character
    if ch == 10
       continue
    else
        fprintf('%c', ch);               % Display character
    end
end

fclose(fileID);
