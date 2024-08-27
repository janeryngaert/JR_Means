%NY_Fed_1bin_tri.m
%This file reads the forecasts records with all probability placed into a
%single bin and fits a triangual distribution to each one. The procedure
%builds upon that used by the NY Fed, but with the major difference that
%the distribution's mode is set to equal the reported point estimate.

%Read the data. (This file gets created by writeForecasts2xls.do.)
data = xlsread('_1bin.xlsx');

%The file is organized with the point estimates in the first column, bins
%1-10 in the next ten columns, and a count of the number of bins with no
%probability mass in the last column.

%Check prerequesites.
assert(size(data,2)==12); %The input data has 12 columns.
assert(~any(ismissing(data))); %No data should be missing.
assert(all(data(:,end)==9)); %Exactly one bin should be not empty.

% Create a vector with the intervals' lower bounds.
ll = [ 12, 8, 4, 2, 0, -2, -4, -8, -12, -38];
% Create a vector with the intervals' upper bounds.
rr = [38, 12, 8, 4, 2, 0, -2, -4, -8, -12];

% The NY Fed sets the bounds of outer bins are set at -38 and 38. I
% think they choose this bound to make the implied mean forcast for those
% filling only one of the outer bins -25 or 25

% Quantiles calculated from each distribution.
desiredQuantiles = [0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95];

%Describe what we are putting in and what we are getting out for each
%record.

% preallocate for mean, l, r, percentiles, 
Epi = zeros(length(p),1);
flag = zeros(length(p),1);
low = zeros(length(p),1);
high = zeros(length(p),1);
quantiles = zeros(length(data),length(desiredQuantiles));

% The NY Fed's convention places the highest inflation outcomes in
% bin1, the next highest in bin2, etc. Therefore, the bin labels
% increment ``backwards'' relative to the inflation outcomes
% themselves.
% To make this code clearer, we flip this ordering of the columns in
% ll, rr, and data.
ll = fliplr(ll);
rr = fliplr(rr);
prob = fliplr(data(:,2:11));

%Umpack the point estimates.
p = data(:,1);
for i = 1:length(data)
        
    % find the single bin with positive probability.
    f = find(prob(i,:) > 0);
    % left bound
    l = ll(f);
    %right bound
    r = rr(f);
    %point estimate
    pp = p(i);




    % save left and right bounds
    % l and r will be used to indicate bounds of support of distribution
    low(i,1) = l;
    high(i,1) = r;
    
    
    % if point estimate falls outside bin with positive probability replace with
    % with the bound of that bin

    if pp > r
        
        pp = r;
        
    elseif pp < l
        
        pp = l;
        
    end


     
    % go through options - if mode (pp) is at a boundary of the
    % distribution, distribution is a right triangle and the pdf is given
    % as the area under a line from l to r
    % want to solve for slope (m) and intercept (b) such that pdf takes
    % positive values and integrates to 1 over l to r


    % CASE ONE
    % start with case where point estimate is at lower bound
    % three vertices of the triangle (T1) are (l,0), (r,0), and (l,h)

    if pp == l
    
        %slope of line between (l,h) and (r,0)
        m = -h/(r - l);
        %intercept
        %use slope intercept form with slope m and point (r,0)
        b = (-1)*r*m;
        
        % pdf function p(x)
        fun1 = @(x) (m.*(x) + b);
        % function p(x)*x - take integral to solve for mean
        fun2 = @(x) (m.*(x.^2)  + b.*x);
        
        %find JR mean
        Epi(i,1) = integral(@(x)fun2(x), l, r);

        %%% calculate percentiles of individual data
        % for case where point estimate is at lower bound, percentile ``centi'' will
        % form a right triangle (T2) with area (1 - cent(j)) between (r,0)
        % (centi, 0) and (centi, height_cent)
        % both centi and height_cent are unknowns, so we will use both the
        % area of T2 and trigonometric ratios between T1 and T2 to solve
        % for height_centi and then centi
        % Area of T2: (1 - cent(j)) = (0.5)*(r- centi)*height_cent
        % Opposite over adjacent for T1 and T2: h/(r - l) = height_cent/(r - centi) -> height_cent = h*(r-centi)/(r-l)
        % Plug height_cent from above into area of T2 formula:
        % (1 - cent(j)) = (0.5)*(r- centi)*[h*(r-centi)/(r-l)]
        % Rearranging
        % (1 - cent(j)) = (0.5)*h*[(r - centi)^2/(r-l)]
        % Plug in h (height of T1)
        % (1 - cent(j)) = (0.5)*[2/(r-l)]*[((r - centi)^2)/(r-l)]
        % Rearranging
        % (1 - cent(j)) = [(r - centi)^2]/[(r-l)^2]
        % r, l, and cent(j) are known and we are looking for centi or the
        % value of the percentile cent(j)
        % Rearrange the above into a
        % quadratic and then use the quadratic formula to solve for the
        % root that is in the bounds [l, r] for each percentile in
        % cent = [0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95]
        % Replacing centi with x here just for quadratic equation:
        % x^2 + -2rx + [r^2 - (1-cent(j))*(r-l)^2]

        for j = 1:length(desiredQuantiles)
        
            % quadratic coefficients and constants from above equation
            aa = 1;
            bb = -2*r;
            cc = (r^2) - (((r - l)^2)*(1 - desiredQuantiles(1,j)));
            % find roots
            d = sqrt(bb^2 - 4*aa*cc);
            y(1) = ( -bb + d ) / (2*aa);
            y(2) = ( -bb - d ) / (2*aa);
            
            %choose root that is in the support of the distribution
            quantiles(i,j) = y(y>l & y<r);
    
        end
        
    % CASE TWO
    % case where point estimate is at upper bound
    % three vertices of the triangle (T1) are (l,0), (r,0), and (r, h)

    elseif pp == r
    
        %slope of line between (l,0) and (r,h)
        m = h/(r - l);
        %intercept - slope intercept form with slope m and point (l,0)
        b = (-1)*l*m;

        % pdf function p(x)
        fun1 = @(x) (m.*(x)  + b);
        % function p(x)*x - take integral to solve for mean
        fun2 = @(x) (m.*(x.^2)  + b.*x);
        Epi(i,1) = integral(@(x)fun2(x), l, r);
            
    %%% calculate percentiles of individual data
    %for case where point estimate is at upper bound, percentile ``centi'' will
    % form a right triangle (T2) with area cent(j) between (l,0)
    % (centi, 0) and (centi, height_cent)
    % both centi and height_cent are unknowns, so we will use both the
    % area of T2 and trigonometric ratios between T1 and T2 to solve
    % for height_centi and then centi
    % Area of T2: cent(j) = (0.5)*(centi - l)*height_cent
    % Opposite over adjacent for T1 and T2: h/(r - l) = height_cent/(centi - l) -> height_cent = h*(centi - l)/(r-l)
    % Plug height_cent from above into area of T2 formula:
    % cent(j) = (0.5)*(centi - l)*[h*(centi - l)/(r-l)]
    % Rearranging
    % cent(j) = (0.5)*h*[(centi - l)^2/(r-l)]
    % Plug in h (height of T1)
    % cent(j) = (0.5)*[2/(r-l)]*[(centi - l)^2/(r-l)]
    % Rearranging
    % cent(j) = [(centi - l)^2]/[(r-l)^2]
    % r, l, and cent(j) are known and we are looking for centi or the    
    % value of the percentile cent(j)
    % Rearrange the above into a
    % quadratic and then use the quadratic formula to solve for the
    % root that is in the bounds [l, r] for each percentile in
    % cent = [0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95]
    % Replacing centi with x here just for quadratic equation:
    % x^2 + -2lx + [l^2 - cent(j)*((r-l)^2)]

        for j = 1:length(desiredQuantiles)
        
            % quadratic coefficients and constants from above equation
            aa = 1;
            bb = -2*l;
            cc = (l^2) - (((r - l)^2)*(desiredQuantiles(1,j)));
            % find roots
            d = sqrt(bb^2 - 4*aa*cc);
            y(1) = ( -bb + d ) / (2*aa);
            y(2) = ( -bb - d ) / (2*aa);
            
            %choose root that is in the support of the distribution
            quantiles(i,j) = y(y>l & y<r);
    
        end
        
    else
    
        % in the last case the point estimate falls between the two bounds,
        % The pdf now forms a triangle between (l,0) (pp,h), and (r,0)
        % The process is similar but now the triangle consists of the area
        % under two lines - L1 between (l,0) and (pp, h) and L2 betwee
        % (pp,h) and (r,0)

        % slope of L1
        m1 = h/(p(i,1) - l);
        % intercept solved for using slope m1 and point (l,0)
        b1 = (-1)*l*m1;

        %slope of L2
        m2 = -h/(r - p(i,1));
        % intercept solved for using slope m2 and point (r,0)
        b2 = (-1)*r*m2;

        % pdf given by L1 between (l,0) and (pp,h)
        fun1 = @(x) (m1.*(x) + b1);
        % pdf given by L2 between (pp,h) and (r,0)
        fun2 = @(x) (m2.*(x) + b2);
        
        %fun3 and fun4 give p(x)*x for both L1 and L2 respectively
        fun3 = @(x) (m1.*(x.^2)  + b1.*x);
        fun4 = @(x) (m2.*(x.^2)  + b2.*x);
        
        %take integral over appropriate ranges to solve for mean

        Epi(i,1) = integral(@(x)fun1(x), l, pp) + integral(@(x)fun2(x), pp, r);  
        
        
        % The point forecast pp splits the distribution into two right
        % triangles that we can use to calculate the percentiles as in the
        % two cases above
        % Call these T3 with vertices (l,0) (pp, 0) and (pp, h) and
        % T4 with vertices (pp, 0), (pp, h) and (r, 0)
        % Both T3 and T4 have the same height as T1 (h)

        %  Solve for the mass of the distribution between (l,0) and (pp,h)

        pr = h*(pp - l)*(0.5);

        %  Any percentile smaller than pr will be calculated similarly to those in CASE TWO
        %  These will use T3 for the trigonometric ratios
        %  Any percentile larger than pr will be calculated similarly to those in CASE ONE
        %  These will use T4 for the trigonometric ratios

        for j = 1:length(desiredQuantiles)
        
            if desiredQuantiles(1,j) < pr

                % variation on percentiles from CASE TWO percentile ``centi'' will
                % form a right triangle (T2) with area cent(j) between (l,0)
                % (centi, 0) and (centi, height_cent)
                % instead of T1 for trigonopetric ratios, we will use the
                % right triangle T3 between (l,0) (pp, 0) and (pp, h)
                % for height_centi and then centi
                % Area of T2: cent(j) = (0.5)*(centi - l)*height_cent
                % Opposite over adjacent for T3 and T2: h/(pp - l) = height_cent/(centi - l) -> height_cent = h*(centi - l)/(pp-l)
                % Plug height_cent from above into area of T2 formula:
                % cent(j) = (0.5)*(centi - l)*[h*(centi - l)/(pp-l)]
                % Rearranging
                % cent(j) = (0.5)*h*[(centi - l)^2/(pp-l)]
                % Plug in h (height of T3)
                % cent(j) = (0.5)*[2/(r-l)]*[(centi - l)^2/(pp-l)]
                % Rearranging
                % cent(j) = [(centi - l)^2]/[(r-l)(pp-l)]
                % r, l, pp and cent(j) are known and we are looking for centi or the    
                % value of the percentile cent(j)
                % Rearrange the above into a
                % quadratic and then use the quadratic formula to solve for the
                % root that is in the bounds [l, pp] for each percentile in
                % cent = [0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95] < pr
                % Replacing centi with x here just for quadratic equation:
                % x^2 + -2lx + [l^2 - (cent(j)*((r-l)(pp-l)))]

                
                %quadratic formula coefficients
                aa = 1;
                bb = -2*l;
                cc = (l^2) - (((pp-l)*(r-l))*(desiredQuantiles(1,j)));
                
                % solve for roots
                d = sqrt(bb^2 - 4*aa*cc);
                y(1) = ( -bb + d ) / (2*aa);
                y(2) = ( -bb - d ) / (2*aa);
            
                % choose value between l and pp
                quantiles(i,j) = y(y >= l & y <= pp);
            
            elseif desiredQuantiles(1,j) >= pr
            
                % variation on percentiles from CASE ONE percentile ``centi'' will
                % form a right triangle (T2) with area (1 - cent(j)) between
                % (centi, 0) and (centi, height_cent) and (r,0)
                % instead of T1 for trigonopetric ratios, we will use the
                % right triangle T4 between (pp, 0) (pp, h) and (r, 0
                % for height_centi and then centi
                % Area of T2: (1 - cent(j)) = (0.5)*(r - centi)*height_cent
                % Opposite over adjacent for T3 and T4: h/(r - pp) = height_cent/(r - centi) -> height_cent = h*(r - centi)/(r - pp)
                % Plug height_cent from above into area of T2 formula:
                % (1 - cent(j)) = (0.5)*(r - centi)*[h*(r - centi)/(r - pp)]
                % Rearranging
                % (1 - cent(j)) = (0.5)*h*[(r - centi)^2/(r - pp)]
                % Plug in h (height of T4)
                % (1 - cent(j)) = (0.5)*[2/(r - l)]*[(r - centi)^2/(r - pp)]
                % Rearranging
                % (1 - cent(j)) = (r - centi)^2/((r - pp)(r-l))


                % quadratic coefficients
                aa = 1;
                bb = -2*r;
                cc = (r^2) - (((r - pp)*(r-l))*(1 - desiredQuantiles(1,j)));
                
                % solve for roots
                d = sqrt(bb^2 - 4*aa*cc);
                y(1) = ( -bb + d ) / (2*aa);
                y(2) = ( -bb - d ) / (2*aa);
            
                % find roots in range [pp, r]
                quantiles(i,j) = y(y >= pp & y <= r);
                
            end
        end
    
    end

    % records point forecast used for calculation of means
    % this could go as (in this case it just replaces the point forecast
    % with the left or right bound in the case that the point forecast
    % falls outside these bounds
    % In cases with more than one bin this can get more complicated

    p2(i,1) = pp;
    
end    
%end



A1 = [data(:, (1:12)), Epi, quantiles, low, high, p2];

xlswrite('Epi_1year_1bin.xlsx', A1);


%%Jane, can I delete the code below? JRC

A = [low, quantiles, high, p, Epi];


% a check for the percentiles - if 
Flag = zeros(length(data), 9);

for i = 1:length(data)
for j = 2:9
if A(i,j) <= A(i, j-1)
Flag(i,j) = 1;
end
end
end

sum(Flag)


