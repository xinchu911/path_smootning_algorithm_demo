clear;
clc;
%标志位:采样使用坐标轴等距->1,其他值使用曲线长等距
flag_axis_interval=0;
%标志位:不画出结果图->0,否则为其他值
flag_plot_result=1;

pc=[0, -0.4, 0.02, -0.004]; %输入三次多项式曲线系数
%初始的x,y均为0
init_x=0;
init_y=0;
knot_interval=5;
knot_num=6;
anchor_num=knot_num-1;
input_points=zeros(knot_num-1,2);
anchor_points=zeros(knot_num,2);
if flag_axis_interval ~= 1
    %采用符号函数求解积分,得等距原曲线采样点
    syms x x0;

    S_length=knot_interval;
    y=pc(1)+pc(2)*x+pc(3)*x^2+pc(4)*x^3;

    f=sqrt(1+diff(y)*diff(y));
    %曲线积分表达式
    for i=1:anchor_num
        S=int(f,x,[init_x,x0]);
        input_points(i,2)=vpasolve(S==(S_length+knot_interval*(i-1)),x0);
        input_points(i,1)=pc(1)+pc(2)*input_points(i,2)+pc(3)*input_points(i,2)^2+pc(4)*input_points(i,2)^3;
    end
    input_points=[[0,0];input_points];
    plot(input_points(1:knot_num,1),input_points(1:knot_num,2),'.','color','r','MarkerSize',10)
    hold on
    %在采样点的中点处算anchor_point
    S_t=int(f,x,[0,x0]);
    anchor_start_y=vpasolve(S==(S_length/2),x0);
    anchor_start_x=pc(1)+pc(2)*anchor_start_y+pc(3)*anchor_start_y^2+pc(4)*anchor_start_y^3;
    for i=1:anchor_num-1
        S_an=int(f,x,[anchor_start_y,x0]);
        anchor_points(i,2)=vpasolve(S_an==(S_length+knot_interval*(i-1)),x0);
        anchor_points(i,1)=pc(1)+pc(2)*anchor_points(i,2)+pc(3)*anchor_points(i,2)^2+pc(4)*anchor_points(i,2)^3;
    end
    anchor_points=[[anchor_start_x,anchor_start_y];anchor_points];
    plot(anchor_points(1:anchor_num,1),anchor_points(1:anchor_num,2),'.','color','b','MarkerSize',10)
    hold on
else
    %按y轴等距采样与anchor_point
    for i=1:knot_num
        input_points(i,2)=(i-1)*knot_interval;
        input_points(i,1)=pc(1)+pc(2)*input_points(i,2)+pc(3)*input_points(i,2)^2+pc(4)*input_points(i,2)^3;
        if i~=1
            anchor_points(i-1,2)=(i-1.5)*knot_interval;
            anchor_points(i-1,1)=pc(1)+pc(2)*anchor_points(i-1,2)+pc(3)*anchor_points(i-1,2)^2+pc(4)*anchor_points(i-1,2)^3;
        end
    end
    plot(input_points(1:knot_num,1),input_points(1:knot_num,2),'.','color','r','MarkerSize',10)
    hold on
    plot(anchor_points(2:anchor_num,1),anchor_points(2:anchor_num,2),'.','color','b','MarkerSize',10)
    hold on
end
%输出采样点anchor_point至csv文件供C++程序读取
No=[1:knot_num].';
title={'No','input_points_x','input_points_y','anchor_points_x','anchor_points_y'};
result_table=table(No,input_points(1:knot_num,1),input_points(1:knot_num,2),double(anchor_points(1:knot_num,1)),double(anchor_points(1:knot_num,2)),'VariableNames',title);
writetable(result_table,'input_data.csv');

y=input_points(1,2):0.01:input_points(knot_num,2);
plot(polyval([-0.004,0.02,-0.4,0],y),y,'b')
hold on

xlabel('x');
ylabel('y');
% legend('分段点','anchor点','输入三项式曲线')

%read结果多项式系数向量,plot结果
if flag_plot_result~=0
    QPSolution=csvread('QPSolution.csv',1,0);
    p=zeros(12,6);
    step=0.01;
    t=0:step:1;
    bounding_box=0.05;
    %画出二次规划后的拟合各段曲线
    for i=1:6
        p(1:6,i)=flip(QPSolution(1+6*(i-1):6*i,:));
        p(7:12,i)=flip(QPSolution(36+1+6*(i-1):36+6*i,:));
        plot(polyval(p(1:6,i),t),polyval(p(7:12,i),t),'g')
        %画出各anchor_point的bounding box
        if i<=5
            y_bouding_low=(anchor_points(i,2)-bounding_box).*ones(2*bounding_box/step+1,1);
            y_bouding_high=(anchor_points(i,2)+bounding_box).*ones(2*bounding_box/step+1,1);
            x_bouding_low=anchor_points(i,1)-bounding_box.*ones(2*bounding_box/step+1,1);
            x_bouding_high=(anchor_points(i,1)+bounding_box).*ones(2*bounding_box/step+1,1);
            plot(anchor_points(i,1)-bounding_box:step:anchor_points(i,1)+bounding_box,y_bouding_low,':','LineWidth',2,'Color','c')
            plot(anchor_points(i,1)-bounding_box:step:anchor_points(i,1)+bounding_box,y_bouding_high,':','LineWidth',2,'Color','c')
            plot(x_bouding_low,anchor_points(i,2)-bounding_box:step:anchor_points(i,2)+bounding_box,':','LineWidth',2,'Color','c')
            plot(x_bouding_high,anchor_points(i,2)-bounding_box:step:anchor_points(i,2)+bounding_box,':','LineWidth',2,'Color','c')
        end
    end
end

