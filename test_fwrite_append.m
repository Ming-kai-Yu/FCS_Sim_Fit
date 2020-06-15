%% test fwrite

A_dat = rand(10, 2);
B_dat = rand(10, 2)+ 1;

fileID = fopen('test.bin','w');
fwrite(fileID, A_dat,'double');
fclose(fileID);

fileID = fopen('test.bin', 'a');
fwrite(fileID, B_dat, 'double');
fclose(fileID);

fileID = fopen('test.bin');
vec = fread(fileID, 40, 'double');
fclose(fileID);

fileID = fopen('test.bin');
ab = fread(fileID, [10, 4], 'double');
fclose(fileID);