function clip(mic_name, mic_name_2, ref_name)
Srate=16000;
t= 25;

file_id=fopen(mic_name, 'r');
x=fread(file_id, inf, 'int16');
fclose(file_id);

s=regexp(mic_name, '\.', 'split');
outfile=strcat(s(1),'_short.',s(2));
outfile=char(outfile);
xfinal=x(1:t*Srate);
file_id=fopen(outfile,'wb');
fwrite(file_id,xfinal,'int16');
fclose(file_id);

file_id=fopen(mic_name_2, 'r');
x2=fread(file_id, inf, 'int16');
fclose(file_id);

s=regexp(mic_name_2, '\.', 'split');
outfile=strcat(s(1),'_short.',s(2));
outfile=char(outfile);
xfinal=x2(1:t*Srate);
file_id=fopen(outfile,'wb');
fwrite(file_id,xfinal,'int16');
fclose(file_id);

file_id=fopen(ref_name, 'r');
r=fread(file_id, inf, 'int16');
fclose(file_id);

s=regexp(ref_name, '\.', 'split');
outfile=strcat(s(1),'_short.',s(2));
outfile=char(outfile);
xfinal=r(1:t*Srate);
file_id=fopen(outfile,'wb');
fwrite(file_id,xfinal,'int16');
fclose(file_id);
end