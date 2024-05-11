% clear all
% mex cec13_func.cpp -DWINDOWS
clear;
clc;
format long;
fun_nums=30;
D=30;
Xmin=-100;
Xmax=100;
nfe_max=10000*D;
runs=51;
%xbest 最优位置矩阵
xbest(runs,D) = inf;
xbest(:) = inf;
%fbest 最优值矩阵
fbest(fun_nums,runs) = inf;
fbest(:) = inf;
%error 矩阵
f_error(fun_nums,runs) = inf;
f_error(:) = inf;
%f_mean 均值矩阵
f_mean(fun_nums) = inf;
f_mean(:) = inf;
%f_median 中值矩阵
f_median(fun_nums) = inf;
f_median(:) = inf;
%f_std 标准差
f_std(fun_nums) = inf;
f_std(:) = inf;
Time = zeros(fun_nums,runs);

targetbest = [100;200;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;
             2000;2100;2200;2300;2400;2500;2600;2700;2800;2900;3000];

fhd=str2func('cec14_func');
fname = ['ADEDMR_',num2str(D),'D.txt'];
f_out = fopen(fname,'wt');
fMedian = ['ADEDMR_',num2str(D),'D_locMedian.txt'];
f_loc = fopen(fMedian,'wt');
%ftime = ['PalmDETime_',num2str(D),'D.txt'];
%f_tout = fopen(ftime,'a');

for i=1:fun_nums
%for i = 14:14 
    fun_num=i;
    disp(['fid:',num2str(fun_num)]);
    fprintf(f_out,'fid:%d\n',fun_num);
    fprintf(f_loc,'fid:%d\n',fun_num);
    for j=1:runs
        [gbest,gbestval,FES,RecordT]= ADEDMR_func(fhd,D,nfe_max,Xmin,Xmax,fun_num,j);
        xbest(j,:)=gbest;   %D维行向量
        fbest(i,j)=gbestval;%单一数值
        disp(['x[',num2str(gbest),']=',num2str(gbestval-targetbest(i),15)]);
        fprintf(f_out,'x[%s]=%s\n',num2str(gbest),num2str(gbestval-targetbest(i)));
        Time(i,j) = RecordT;
    end
    [bestval,best] = min(fbest(i,:));
    loc = xbest(best,:);
    disp(['Best[',num2str(loc),']=',num2str(bestval-targetbest(i),15)]);
    fprintf(f_out,'Best[%s]=%s\n',num2str(loc),num2str(bestval-targetbest(i)));
    %distortion, accuracy loss
    %f_mean(i)=mean(fbest(i,:));
    %disp(['mean[ ',num2str(i),']=',num2str(f_mean(i)-targetbest(i),15)]);
    f_error(i,:) = fbest(i,:) - targetbest(i);
    f_mean(i)=mean(f_error(i,:));
    f_median(i) = median(f_error(i,:));
    [value,index] = sort(f_error(i,:));
    if mod(runs,2)==0
        locMedian = [index(runs/2),index(runs/2+1)];
    else
        locMedian = index(ceil(runs/2));
    end
    fprintf(f_loc,'LocMedian = %s\n',num2str(locMedian));
    f_std(i) = std(f_error(i,:));
    disp(['mean[',num2str(i),']=',num2str(f_mean(i),15)]);
    fprintf(f_out,'mean[%s]=%s\n',num2str(i),num2str(f_mean(i)));
    disp(['median[',num2str(i),']=',num2str(f_median(i),15)]);
    fprintf(f_out,'median[%s]=%s\n',num2str(i),num2str(f_median(i)));
    disp(['std[',num2str(i),']=',num2str(f_std(i),15)]);
    fprintf(f_out,'std[%s]=%s\n',num2str(i),num2str(f_std(i)));
    %{
    for j = 1: runs
    	disp(['T2=',num2str(Time(i,j),15)]);
    	fprintf(f_tout,'T2=\t%.15f\n',Time(i,j));
    end;
    MeanT=mean(Time(i,:));
    disp(['meanT[',num2str(i),']=',num2str(MeanT,15)]);
    fprintf(f_tout,'MeanT2=\t%.15f\n',MeanT);
    %}
end
fclose(f_out);
fclose(f_loc);
%fclose(f_tout);
clear all;