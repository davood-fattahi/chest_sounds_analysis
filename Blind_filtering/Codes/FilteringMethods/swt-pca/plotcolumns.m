function plotcolumns(varargin)

if length(varargin)==1
    Y=varargin{1};
    J=size(Y,2);
    for j=1:J
            subplot(J,1,j)
            plot(Y(:,j));
            grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end
elseif length(varargin)>=2 && ismatrix(varargin{2})
    Y=varargin{1};
    X=varargin{2};
    J=size(Y,2);
    for j=1:J
            subplot(J,1,j)
            plot(X,Y(:,j),varargin{3:end});
            grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end
elseif length(varargin)>=2 && ~ismatrix(varargin{2})
    Y=varargin{1};
    J=size(Y,2);
    for j=1:J
            subplot(J,1,j)
            plot(Y(:,j),varargin{2:end});
            grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end 
end


