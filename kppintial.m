% k++ initialization outlined in project instructions (doing what kmeans does without using kmeans function in MATLAB).


load Q1data.mat;

function [c,IndexSet]= kppinitial(XData,k,varargin)

[n,m] = size(XData);

IndexSet = randi(k,n,1);

c = zeros(k,m);

if nargin == 2
    randIndex = randi(n);

    c(1,:) = XData(randIndex,:);
else
    init = varargin{1};
    c(1,:) = XData(init,:);
end

closestCluster=zeros(n,1);

for l = 2:k 
    for d=1:n

        % Store the coordinates of the current data vector
        xD = XData(d,:);

        % Set the minimum distance tracker to be a very large number
        sqDistMin=1e16;

        % Find the closest weight vector (cluster) to the current data
        % vector
        for i=1:l-1
            sqDist = norm(c(i,:)-xD,2);

            % If the distance is less than the current min, assign the
            % current data vector to this cluster
            if sqDist<sqDistMin
                closestCluster(d)=i;
                sqDistMin=sqDist;
            end

        end
    end

     % Update the assignments of the data vectors to their new clusters
    IndexSet = closestCluster;

    DistClust = zeros(l-1,1+m);

    for y = 1:l-1
        PntClosest = XData(IndexSet == y,:);

        [MaxDistFromClosestCluster,IndexClosest] = max(sum(PntClosest - c(y,:)).^2,2);
        DistClust(y,:) = [MaxDistFromclosestCluster,PntClosest(IndexClosest,:)];
    end

    [~, NextCentroid] = max(DistClust(:,1));
    c(l,:) = DistClust(NextCentroid,2:end);
end

    for d=1:n

        % Store the coordinates of the current data vector
        xD = XData(d,:);

        % Set the minimum distance tracker to be a very large number
        sqDistMin=1e16;

        % Find the closest weight vector (cluster) to the current data
        % vector
        for i=1:k
            sqDist = norm(c(i,:)-xD,2);

            % If the distance is less than the current min, assign the
            % current data vector to this cluster
            if sqDist<sqDistMin
                closestCluster(d)=i;
                sqDistMin=sqDist;
            end

        end
    end

    % Update the assignments of the data vectors to their new clusters
    IndexSet = closestCluster;
end
