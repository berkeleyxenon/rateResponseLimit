function [x,y] = makeSimple2Dcontour(tmpx,tmpy)
% little function to sort vector for making contour plot...
% X,Y are a vectors of scatter plot points
% the goal is to order them so that plot(x,y,'b-') returns a reasonable result
	x = tmpx(1); tmpx(1)=-1;
	y = tmpy(1); tmpy(1)=-1;
	for ii=1:length(tmpx)-1 
		de=sqrt( (tmpx-x(end)).^2 + (tmpy-y(end)).^2); % euclidean distance
		nI = find( de==(min(de(de>0))) ); 
%plot(x,y,'go-');
 %       if ii==38;keyboard;end;

		x = [x tmpx(nI)]; tmpx(nI) = -1;
		y = [y tmpy(nI)]; tmpy(nI) = -1;
	end

%x(end+1) = x(1);
%y(end+1) = y(1);
