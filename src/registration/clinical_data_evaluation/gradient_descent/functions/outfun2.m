function stop = outfun2(x, optimValues, state) % output function

global save_folder

stop = false;

fileName = strcat(save_folder,'/optim_output.txt');
[fid, message] = fopen(fileName,'a');

switch state
    
    case 'init'
        
        if fid==-1
            error('globaloptim:psoutputfile:fileError','Error trying to write to %s:\n%s',fileName,message)
        end
        
    case 'iter'
        
        Iter = optimValues.iteration;
        Fcount   = optimValues.funccount;
        Fval     = optimValues.fval;
        Step     = optimValues.stepsize;
        Firstopt = optimValues.firstorderopt;
        % Write to the file
        fprintf(fid,'Iter: %d \t Func-count: %d \t f(x): %5.2f \t Step-size: %5.5f \t First-order-opt: %5.2f \n',Iter, Fcount, Fval, Step, Firstopt);

    case 'done'
        
        Iter     = optimValues.iteration;
        Fcount   = optimValues.funccount;
        Fval     = optimValues.fval;
        Step     = optimValues.stepsize;
        Firstopt = optimValues.firstorderopt;
        % Write to the file
        fprintf(fid,'Iter: %d \t Func-count: %d \t f(x): %5.2f \t Step-size: %5.5f \t First-order-opt: %5.2f \n',Iter, Fcount, Fval, Step, Firstopt);

        fprintf(fid,'\nOptimization terminated.\n');
        fclose(fid);
        
end

