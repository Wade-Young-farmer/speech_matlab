function clip(file_name)
s=regexp(file_name, '\.', 'split');
outfile=strcat(s(1),'_short.',s(2));
outfile=char(outfile);
if strcmp(s(2),'pcm')
    Srate=16000;
    t= 25;

    file_id=fopen(file_name, 'r');
    x=fread(file_id, inf, 'int16');
    fclose(file_id);
    
    xfinal=x(1:t*Srate);
    file_id=fopen(outfile,'wb');
    fwrite(file_id,xfinal,'int16');
    fclose(file_id);
elseif strcmp(s(2), 'wav')
    [x, fs] = audioread(file_name);
    x=x(fs*9 + 1:fs*14,:);
    audiowrite(outfile,x,fs);
end
end