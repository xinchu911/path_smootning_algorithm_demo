QPSolution=[     -1.4153
1.55627
-0.949615
0.00933605
-0.078476
0.0175976
-0.860235
-0.540929
-1.21643
-0.128592
0.145928
-0.0174127
-2.61768
-2.86291
-0.900703
0.280994
-0.205972
0.000483498
-6.30581
-4.64275
-1.28867
-0.53806
0.0499772
0.0230293
-12.7023
-8.51917
-2.37268
-0.107859
0.0469612
-0.0108594
-23.6659
-13.4546
-2.52308
-0.028608
0.0093396
-0.000972268
3.59892
-3.89078
2.75151
-0.0119324
0.0133618
-0.0247237
2.43639
1.50616
2.54868
-0.205722
-0.360552
0.0501808
5.97505
4.79501
0.270095
-1.14612
0.438078
0.00232172
10.3344
3.66082
-0.516519
0.629408
-0.100232
-0.0443678
13.9635
3.89328
0.326676
-0.215198
-0.0564612
0.0271413
 17.939
3.81094
-0.386252
-0.169629
0.0366777
0.00254652];
p=zeros(12,6);
t=0:0.01:1;
bounding_box=0.05;
input_points_x=[0, -1.83034258887661, -4.92346969167710, -10.1590713469645, -19.7044036200127, -36.8875134535896];
input_points_y =[0, 4.65129686919330, 9.02708011754813, 12.7606580112594, 16.6910146745120, 21.0541698358764];
verification_points_x =[-0.877987143447578, -1.44278900367223, -4.25332620694243, -9.01559868875152, -17.5712730896452];
verification_points_y =[2.340665195059720, 3.78122865939735, 8.32673718677687, 12.1104852113046, 15.9671921248409];
for i=1:6
%     p(1:6,i:i)=flip(QPSolution(1+6*(i-1):6*i,:));
%     p(7:12,i:i)=flip(QPSolution(36+1+6*(i-1):36+6*i,:));
%     plot(polyval(p(1:6,i:i),t),polyval(p(7:12,i:i),t),'g')
    %knot分段点
    plot(input_points_x(i),input_points_y(i),'.','color','r','MarkerSize',10)
    %anchor points
    if i<=5
        plot(verification_points_x(i),verification_points_y(i),'.','color','b','MarkerSize',10)
    end
%     if i~=1
        %plot(polyval(p(1:6,i:i),0.5),polyval(p(7:12,i:i),0.5),'.','color','r','MarkerSize',10)
%         plot(polyval([-0.004,0.02,-0.4,0],polyval(p(7:12,i:i),0.5)),polyval(p(7:12,i:i),0.5),'.','color','m','MarkerSize',10)
%     end
%     if i<=5
%             y_bouding_low=(verification_points_y(i)-bounding_box).*ones(2*bounding_box/0.01+1,1);
%             y_bouding_high=(verification_points_y(i)+bounding_box).*ones(2*bounding_box/0.01+1,1);
%             x_bouding_low=verification_points_x(i)-bounding_box.*ones(2*bounding_box/0.01+1,1);
%             x_bouding_high=(verification_points_x(i)+bounding_box).*ones(2*bounding_box/0.01+1,1);
%             plot(verification_points_x(i)-bounding_box:0.01:verification_points_x(i)+bounding_box,y_bouding_low,':','LineWidth',2,'Color','c')
%             plot(verification_points_x(i)-bounding_box:0.01:verification_points_x(i)+bounding_box,y_bouding_high,':','LineWidth',2,'Color','c')
%             plot(x_bouding_low,verification_points_y(i)-bounding_box:0.01:verification_points_y(i)+bounding_box,':','LineWidth',2,'Color','c')
%             plot(x_bouding_high,verification_points_y(i)-bounding_box:0.01:verification_points_y(i)+bounding_box,':','LineWidth',2,'Color','c')
%     end
    hold on
end
%输入三次曲线
y=0:0.01:21.1;
plot(polyval([-0.004,0.02,-0.4,0],y),y,'b')
