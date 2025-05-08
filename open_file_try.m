clear;
clc;
filename = 'bablo.txt';
try
    fileID = fopen("bablo.txt",'w');
    if fileID == -1
        error('FileDoesNotExist:FileOpenError','File open error')
    end

catch ME
    disp(ME.identifier)
end