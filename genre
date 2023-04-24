function [ge]=genre(x)   %Code to be stored in a file named 'genre.m'
for i=1:(length(x)-1)
    for j=(i+1):length(x)
	    x(j,:)=x(j,:)-((x(j,i)/x(i,i))*x(i,:));
	end
end
ge=x;