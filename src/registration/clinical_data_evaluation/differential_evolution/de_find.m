 function V = de_find(SwarmR, target, outer_shift, min_bounds, max_bounds, best_template, pen)
    persistent img;
    persistent outer_shift_cart;
    persistent template;
    persistent lb;
    persistent ub;
    persistent penalties;
    if nargin > 1
        outer_shift_cart = outer_shift;
        img = target;%imrotate(target,-90);
        template = best_template;
        lb = min_bounds;
        ub = max_bounds;
        penalties = pen;
    else
       [SwarmSize, d1] = size(SwarmR);
%        n_points = d1/2;
%        innerSwarm = zeros(n_points, 2);
%        movements_cart = zeros(size(innerSwarm));
       appD = zeros(SwarmSize,1);
       for q=1:SwarmSize
           appD(q) = main_mex(img,SwarmR(q,:),lb,ub,template,ceil(outer_shift_cart),penalties);
%            inner_move = reshape(SwarmR(q,:),n_points,2).*(ub-lb)+lb;
%            [movements_cart(:,1) movements_cart(:,2)] = pol2cart(inner_move(:,1),inner_move(:,2));
%            innerSwarm(1,:) = ceil(movements_cart(1,:));
%            for i=2:size(innerSwarm,1)
%                innerSwarm(i,:) = ceil(innerSwarm(i-1,:) + movements_cart(i,:));
%            end
%            innerSwarm(innerSwarm<1) = 1;innerSwarm(innerSwarm>size(img,1)) = size(img,1);
%            outerSwarm = ceil(innerSwarm+outer_shift_cart);
%            outerSwarm(outerSwarm<1) = 1;outerSwarm(outerSwarm>size(img,1)) = size(img,1);
%            appD(q) = pso_hippo_fitness([innerSwarm; outerSwarm],img,template,penalties);
       end
       V = -appD; % minimize
    end
 end

