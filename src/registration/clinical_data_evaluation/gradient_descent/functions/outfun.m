function stop = outfun(x, optimValues, state) % output function

global save_folder

stop = false;
grad = optimValues.gradient;

fileName = strcat(save_folder,'/grad_output.txt');
[fid, message] = fopen(fileName,'a');

switch state
    
    case 'init'
        
        if fid==-1
            error('globaloptim:psoutputfile:fileError','Error trying to write to %s:\n%s',fileName,message)
        end
        
    case 'iter'
        
        Iter = optimValues.iteration;
        % Write to the file
        fprintf(fid,'Iter: %d \t',Iter);
        for i=1:length(grad)
            fprintf(fid,'%5.4f \t', grad(i));
        end
        fprintf(fid,'\n');
        fprintf(fid,'%5.4f \n', norm(grad));

    case 'done'
        
        Iter = optimValues.iteration;
        % Write to the file
        fprintf(fid,'Iter: %d \t',Iter);
        for i=1:length(grad)
            fprintf(fid,'%5.4f \t', grad(i));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\nOptimization terminated.\n');
        fclose(fid);
        
end

