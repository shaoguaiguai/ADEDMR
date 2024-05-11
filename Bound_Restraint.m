%This function is used for LS bound checking 
function p = Bound_Restraint(p,lowB,upB)
    for i = 1 : size(p,1)
        upper = double(gt(p(i,:),upB));
        lower = double(lt(p(i,:),lowB));
        up = find(upper == 1);
        lo = find(lower == 1);
        if (size(up,2)+ size(lo,2) > 0 )
            for j = 1 : size(up,2)
%                 fprintf('here');
                p(i, up(j)) = rand*(upB(up(j)) - lowB(up(j)))+ lowB(up(j));
            end
            for j = 1 : size(lo,2)
                p(i, lo(j)) = rand*(upB(lo(j)) - lowB(lo(j)))+ lowB(lo(j));
            end
        end
    end
end