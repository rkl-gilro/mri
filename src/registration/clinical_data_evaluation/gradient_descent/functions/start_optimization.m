function start_optimization(  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Start optimization calculates the barycentric coordinates of every 
%%  intersection points given the initial triangulation 'source_tri'
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global optimizer
global t_reg

global target_tri_ax
global target_tri_sag
global target_tri_cor

global source_tri

if optimizer.metaheuristic
    
    if optimizer.views == 2
        
        name_func = 'myfun_unc2';
        %% Call the optimization function
        tic
        [xfinal_v, fval] = optimizationWithDE(1, 6*size(source_tri.X,1),[],[],[], [], [], [], [], [], name_func);
        t_reg = toc
        
        xfinal = reshape(xfinal_v', 2*size(source_tri.X,1),3) +  [source_tri.X;source_tri.X];
        
        target_tri_ax  = TriRep(source_tri.Triangulation, xfinal(1:size(source_tri.X,1),:));
        target_tri_sag = TriRep(source_tri.Triangulation, xfinal(size(source_tri.X,1)+1:end,:));
        
    end
    if optimizer.views == 3
        
        disp('--------- Differential Evolution  -----')
        disp(' ')
        name_func = 'myfun_meta';
        %% Call the optimization function
        tic
        [xfinal, fval] = optimizationWithDE(1, 6*size(source_tri.X,1),[],[],[], [], [], [], [], [], name_func);
        t_reg = toc
        
        mesh1 = xfinal(1:length(xfinal)/3);
        mesh2 = xfinal(length(xfinal)/3 + 1 :2 * length(xfinal)/3 );
        mesh3 = xfinal(2*length(xfinal)/3 + 1 :end );
        
        mesh11 = reshape(mesh1',2,length(xfinal)/6);
        mesh21 = reshape(mesh2',2,length(xfinal)/6);
        mesh31 = reshape(mesh3',2,length(xfinal)/6);
        
        %% Define the triangulations
        target_tri_ax  = TriRep(source_tri.Triangulation, [mesh11' + source_tri.X(:,1:2) source_tri.X(:,3)]);
        target_tri_sag = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh21' + source_tri.X(:,2:3)]);
        target_tri_cor = TriRep(source_tri.Triangulation, [mesh31(1,:)' + source_tri.X(:,1) source_tri.X(:,2) mesh31(2,:)' + source_tri.X(:,3)]);
        
    end
    
    
end

if optimizer.steepestdescent
    
    if optimizer.views == 2
        
        [in, ~] = preparing4;
        
        %% Options
        if optimizer.grad
            opt4 = optimset('Display','iter','Jacobian','on','DerivativeCheck',optimizer.checkderiv,'MaxIter', optimizer.maxiter);
        else
            opt4 = optimset('Display','iter','MaxIter', optimizer.maxiter);
        end
        
        %% Define the initial condition and the boundaries
        tmp = source_tri.X';
        tmp2 = tmp(:);
        X0 = [tmp2' tmp2'];
        lb = -20.*ones(size(X0)) + X0;
        ub =  20.*ones(size(X0)) + X0;
        
        %% Call the optimization function
        tic
        [xfinal, ~ , exitflag, ~] = lsqnonlin(@(t)myfun4(t,in),X0', lb', ub', opt4);
        t_reg = toc
        
        mesh1 = xfinal(1:length(X0)/2);
        mesh2 = xfinal(length(X0)/2 + 1 :end);
        
        mesh11 = reshape(mesh1',3,length(X0)/6);
        mesh21 = reshape(mesh2',3,length(X0)/6);
        
        target_tri_ax  = TriRep(source_tri.Triangulation, [mesh11' source_tri.X(:,3)]);
        target_tri_sag = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh21']);
        
    end
    
    if optimizer.views == 3
        
        [~, ~] = preparing4;
        
        if optimizer.grad
            disp('--------- Gradient Descent, gradient provided -----')
            options = optimset('LargeScale','off','GradObj','on','TolFun', optimizer.tolfun,'TolX', optimizer.tolx,'Display','iter','MaxIter', optimizer.maxiter,'OutputFcn', {@outfun,@outfun2});
        else
            disp('--------- Gradient Descent, no gradient provided -----')
            options = optimset('LargeScale','off','TolFun', optimizer.tolfun,'TolX', optimizer.tolx, 'Display','iter','MaxIter', optimizer.maxiter,'OutputFcn', {@outfun,@outfun2});
        end
        
        %% Define the initial condition
        tmp1  = source_tri.X(:,1:2)';
        tmp12 = tmp1(:);
        tmp2  = source_tri.X(:,2:3)';
        tmp22 = tmp2(:);
        tmp3  = [source_tri.X(:,1) source_tri.X(:,3)]';
        tmp32 = tmp3(:);
        X0 = [tmp12' tmp22' tmp32'];
        
        lb = -20.*ones(size(X0)) + X0;
        ub =  20.*ones(size(X0)) + X0;
        
        %% Call the optimization function
        tic
        if ~optimizer.con
            [xfinal, ~ , exitflag, ~] = fminunc(@myfun_unc_orthoJ, X0, options);
        else
            [xfinal, ~ , exitflag, ~] = fmincon(@myfun_unc_orthoJ, X0, [],[],[],[], lb, ub, [], options);
        end
        t_reg = toc
        
        exitflag
        mesh1 = xfinal(1:length(xfinal)/3);
        mesh2 = xfinal(length(xfinal)/3 + 1 :2 * length(xfinal)/3 );
        mesh3 = xfinal(2*length(xfinal)/3 + 1 :end );
        
        mesh11 = reshape(mesh1',2,length(xfinal)/6);
        mesh21 = reshape(mesh2',2,length(xfinal)/6);
        mesh31 = reshape(mesh3',2,length(xfinal)/6);
        
        %% Define the triangulations
        target_tri_ax  = TriRep(source_tri.Triangulation, [mesh11' source_tri.X(:,3)]);
        target_tri_sag = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh21']);
        target_tri_cor = TriRep(source_tri.Triangulation, [mesh31(1,:)' source_tri.X(:,2) mesh31(2,:)']);
        
    end
end


% disp('--------- Plot the data points  -----')

%% Plot the data points  
% plot3(var_array1(:,1),var_array1(:,2),var_array1(:,3),'r*','MarkerSize',2);hold on
% plot3(var_array2(:,1),var_array2(:,2),var_array2(:,3),'b*','MarkerSize',2);hold on
% plot3(var_array3(:,1),var_array3(:,2),var_array3(:,3),'k*','MarkerSize',2);hold on


%% Clear memory

varlist = {'X_ax','Y_ax','Z_ax','X_sag','Y_sag','Z_sag','M','rows','cols','total_ax','total_sag','N1','N2','ax','sag','x','y','z','options','i','j','k','A1','b1','Aeq1','beq1','A2','b2','Aeq2','beq2',...
           'A','b','Aeq','beq','x0','lambda','V1','nr','ne','vd','ind_tmp','ind_tmp2','ind_tmp3','ind_tmp4','ind_tmp5','ind_tmp6','ind_tmp7','ind_tmp8','tmp_v1','tmp_v',...
           'new_i','new_j','new_k','neig','A_c','b_c','Aeq_c','beq_c','bb','nx','ny','nz','tmp','s2ind','l_x','l_y','l_z','vari','dt','trep'};
clear(varlist{:})
clear varlist


