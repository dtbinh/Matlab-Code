% This script will edit the bibTeX file to remove the useless and
% troublesome fields

fid_in = fopen('../../Dropbox/Literature/Writing/THESIS/AMS LaTeX Package v4.3.1/refs.bib','r+');
fid_out = fopen('../../Dropbox/Literature/Writing/THESIS/AMS LaTeX Package v4.3.1/refs_edited.bib','w');

while (~feof(fid_in))
    s = fgets(fid_in);
    a=strfind(s,'url = ');
    b=strfind(s,'doi = ');
    c=strfind(s,'note = ');
    string_found=length([a,b,c]);
    if string_found
        % Do Nothing
    else
        fprintf(fid_out,'%s',s);
    end
end


fclose('all');