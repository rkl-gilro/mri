function [x,y,z] = calculate4corners( M, i , j )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculates the four corners of the plane determine by the matrix M
%% From      pixel --> RCS
%%
%% Inputs:  1. M -> transform matrix (4 x 4)  ( pixel --> RCS )
%% Outputs: 1. x -> array of x coordinates of the four points 
%%          2. y -> array of y coordinates of the four points
%%          3. z -> array of z coordinates of the four points
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i = [0 511];
% j = [0 511];

x = zeros(1,length(i)*length(j));
y = zeros(1,length(i)*length(j));
z = zeros(1,length(i)*length(j));


for k=1:length(j)
    for l=1:length(i)
        p = M * [j(k) i(l) 1]'; % 0
        
        x((k-1)*length(j) + l) = p(1);
        y((k-1)*length(j) + l) = p(2);
        z((k-1)*length(j) + l) = p(3);
        
    end
end
P1 = [x(1) y(1) z(1)];
P2 = [x(2) y(2) z(2)];
P3 = [x(3) y(3) z(3)];
P4 = [x(4) y(4) z(4)];

x = [P3(1) P4(1) P2(1) P1(1)];  
y = [P3(2) P4(2) P2(2) P1(2)];   
z = [P3(3) P4(3) P2(3) P1(3)];

