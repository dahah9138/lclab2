function fread_foel_to_mat(file)
    [p,f,~]=fileparts(file);
    data = fread_foel(file);
    nn = data{10,2};
    fnew = [p,'/',f,'.mat'];
   
    save(fnew,'nn');
end