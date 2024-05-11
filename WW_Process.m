%This function is used as a local search, and creates some
%new points based on Wavelet Walk.


%The input function is:                                                   %
%Point: the input point which is going to be diffused                     %
%S: structure of problem information                                      %
%g: generation number                                                     %
%BestPoint: the best point in group                                       %                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The output function is:                                                  %
%createPoint: the new points created by Diffusion process                 %
%fitness: the value of fitness function                                   %


function [createPoint, fitness] = WW_Process(Point,s,g,BestPoint,lu,fhd)
     if g==1
          GeneratePoint = Point+exp(1)/sqrt(5).*pi^(-1/4).*((Point-BestPoint).^2).*exp(-(Point-BestPoint).^2/2);

     else
%          GeneratePoint = BestPoint + 0.5*(abs((Point-BestPoint))) * tan(pi * (rand - 0.5))+(randn*BestPoint - randn*Point) ;
%          GeneratePoint = cos(nfes*pi/3)*BestPoint+(randn*BestPoint - randn*Point);
       GeneratePoint =   exp(1)/sqrt(5).*pi^(-1/4).*((Point-BestPoint).^2).*exp(-(-Point+BestPoint).^2/2);
     end
% 
%     else
%        pra=sqrt(1+nfes)/(nfes);
%        GeneratePoint = normrnd(BestPoint, pra*(abs((Point-BestPoint))), [1 size(Point,2)])+(randn*BestPoint - randn*Point) ;
%     end
     
    
    
    %check bounds of generated point
     GeneratePoint = Bound_Restraint(GeneratePoint,lu(1,:),lu(2,:));
    
%     size(GeneratePoint)  
%     for i=1:size(Point,2)
%         if GeneratePoint(1,i) > S.Uband
%             fprintf('violate upper');
%         end
%         if GeneratePoint(1,i) < S.Lband
%              fprintf('violate lower');
%         end  
%     end
    
%     fitness = cec14_func(GeneratePoint',s); 
    fitness = feval(fhd,GeneratePoint',s);

    createPoint = GeneratePoint;
    %======================================================================
end