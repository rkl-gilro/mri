lista=dir('*.c');
for i=1:length(lista)    
    mex (lista(i).name);
end
