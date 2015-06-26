function y= draw_dec_surf_svm(XS, NORMAL, DTA, kernel_function)
% Figura con el hiperplano de decisión.
P = length(NORMAL);%+length(DTA);
labels= ones(P,1);  labels(NORMAL)= -1;

plot3(XS(NORMAL,1),XS(NORMAL,2),XS(NORMAL,3),'b*'); hold on;
h= plot3(XS(DTA,1),XS(DTA,2),XS(DTA,3),'rs'); hold off;
grid;
hAxis = get(h,'parent');
lims = axis(hAxis);

N=100;
[Xc, Yc, Zc]= meshgrid(linspace(lims(1),lims(2),N),linspace(lims(3),lims(4),N),linspace(lims(5),lims(6),N));

output= zeros(N,N,N);

tr_data= XS(:,[ 1 2 3 ]);
train= 1:P;
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','linear');
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','quadratic');
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','polynomial');
svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function',kernel_function);

for xc=1:N
    for yc=1:N
        for zc=1:N
            output(xc,yc,zc)= eval_svmStruct(svmStruct,[Xc(xc,yc,zc) Yc(xc,yc,zc) Zc(xc,yc,zc)]);
        end
    end
end

plot3(XS(NORMAL,1),XS(NORMAL,2),XS(NORMAL,3),'bo','MarkerSize',7); hold on;
h= plot3(XS(DTA,1),XS(DTA,2),XS(DTA,3),'rs','MarkerSize',7);
hAxis = get(h,'parent');
lims = axis(hAxis);

% Representamos los vectores de soporte.
hold on;
sv = svmStruct.SupportVectors;
% see if we need to unscale the data
scaleData = svmStruct.ScaleData;
if ~isempty(scaleData)
    for c = 1:size(sv, 2)
        sv(:,c) = (sv(:,c)./scaleData.scaleFactor(c)) - scaleData.shift(c);
    end
end
% plot support vectors
hSV = plot3(sv(:,1),sv(:,2),sv(:,3),'*k');

hpatch = patch(isosurface(Xc,Yc,Zc,output,0));
isonormals(Xc,Yc,Zc,output,hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
view([-126,36]);
axis tight
camlight left; 
set(gcf,'Renderer','zbuffer'); 
lighting phong;
grid;
hold off;

legend('NOR','AD','Support vectors','Decision surface');
title(['#SVs= ' num2str(size(sv,1))]);
y= 1;
