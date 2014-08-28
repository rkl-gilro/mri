function [INT,Convert]=fac2int(R)
% function [INT,Convert]=fac2int(R)
% fac2int deals with string entries in data matrices 
% for plots and pivottables
% R=either Matrix,Cell array of string elements, or Cell array of mixed elements
% INT: Integer representation of factor 
% Convert.isnum = [vector of falgs if the column was numeric]
% Convert.names{} = if string, it stores the names for the codes.
if isnumeric(R) | islogical(R)
    INT=R;
    Convert.isnum=ones(1,size(R,2));
elseif iscell(R)
    [rows,cols]=size(R);
    if (rows>1)
        R={R};
        [rows,cols]=size(R);
    end;
    INT=zeros(size(R{1},1),cols);
    for i=1:cols
        if iscell(R{1,i});
            Convert.isnum(i)=0;
            [Convert.names{i},d,INT(:,i)]=unique(R{1,i});
        else
            Convert.isnum(i)=1;
            INT(:,i)=R{1,i};
        end;
    end;
end;