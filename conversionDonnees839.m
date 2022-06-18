clear;clc;close all;
tic
fid = fopen('D:\Backup Dell\Dell-D\Enseignement\Encadrement\2021-2022\ASAL\839\R-REC-P.839-4-201309-I!!ZIP-E\h0.txt');
h0_txt = fread(fid);
fclose(fid);

k=1;
value=[];
l=1;
c=1;
for i=1:length(h0_txt)
    if double(h0_txt(i)) == 9
        value = str2num(value);
        h0(l,c)=value;
        c=c+1;
        value=[];
    elseif double(h0_txt(i)) == 13
        value = str2num(value);
        h0(l,c)=value;
        value = [];
        i=i+1;
        l=l+1;
        c=1;
    elseif double(h0_txt(i)) == 10
    else
        value = [value,char(h0_txt(i))];
    end
end

save('h0.mat','h0')