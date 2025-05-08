clc; clear;

fileID = -1;  

try
    fileID = fopen('tsunamis.txt');
    if fileID == -1
        error('FileOpenError:FileNotFound', 'you cant modify this file');
    end
    tline = fgetl(fileID);
    fclose(fileID);
catch ME
    % disp("Something went wrong:");
    disp(ME.identifier)
    % disp(ME.message);
end

% clc; clear;
% 
% filename = 'output1.txt';  % Name of the new file
% 
% % Open file for writing ('w' mode creates a new file or overwrites an existing one)
% fileID = fopen(filename, 'w+');
% 
% if fileID == -1
%     error('FileError:CannotOpenFile', 'Failed to open the file for writing.');
% end
% 
% % Write some content to the file
% 
% fprintf(fileID, 'This is a test.\n');
% fprintf(fileID, 'Writing numbers: %d, %f\n', 42, 3.14159);
% 
% % Always close the file when done
% fclose(fileID);
% 
% disp('File writing completed successfully.');
