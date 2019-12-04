function clip(pcm_name)
Srate=16000;
t= 25;

file_id=fopen(pcm_name, 'r');
x=fread(file_id, inf, 'int16');
fclose(file_id);

s=regexp(pcm_name, '\.', 'split');
outfile=strcat(s(1),'_short.',s(2));
outfile=char(outfile);
xfinal=x(1:t*Srate);
file_id=fopen(outfile,'wb');
fwrite(file_id,xfinal,'int16');
fclose(file_id);
end