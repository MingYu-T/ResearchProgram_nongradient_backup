
function Phi_s = generateTDF(Bx,By,Br)
% - Global Variables
    global numCPoint numSpline ord DH DW nelx nely 

% - Parameters of Splines
    n = numCPoint-1;                    % n+1 is the number of control points
    k = ord;                              % order of the polynomial, i.e. k=degree+1
    [t,Range] = KnotVector(k,n); % knot vector
    SampleSize = 500;
    % Sample Size is the number of divisions between the knot values
    Px = zeros(numSpline,SampleSize*(n-1));
    Py = zeros(numSpline,SampleSize*(n-1));
    r =Br;
    

% - Parameters of Level Set Function 
    [X,Y] = meshgrid(0:DW/nelx:DW,0:DH/nely:DH);
    [ysize,xsize] = size(X);
    Phi = zeros(ysize,xsize,numSpline);

% - Base Function
    x = zeros(1,(n+k)*SampleSize+1);
    x(((k-1)*SampleSize+1):((n+1)*SampleSize+1)) = 0:1/SampleSize:Range;
    x(((n+1)*SampleSize+2):end) = Range;
    x = x./Range; % normalize

    N1 = FirstOrderBSplineFunction(k,t,x,SampleSize);
    N = GeneralOrderBsplineFunction(k,t,x,SampleSize,N1);

% - B-spline points
    for i = 1:numSpline
        Px(i,:) = Bx(i,:)*N(1:n+1,(k-1)*SampleSize+1:(n+1)*SampleSize);
        Py(i,:) = By(i,:)*N(1:n+1,(k-1)*SampleSize+1:(n+1)*SampleSize);
    end

% - Level Set Function
    for s = 1:numSpline
        for i = 1:ysize
            for j = 1:xsize
                dmin = 100;
                for k = 1:SampleSize*(n-1)
                    d = sqrt((Px(s,k)-X(i,j))^2+(Py(s,k)-Y(i,j))^2);
                    if d < dmin
                        dmin = d;
                    end
                end
                Phi(i,j,s) = r(s)-dmin;   % Phi of the splines
            end
        end
    end
    Phi_s = Phi(:,:,1);
    for s = 2:numSpline
        Phi_s = max(Phi_s,Phi(:,:,s));   % Phi of the whole structure
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knot Vector
function [t,Range] = KnotVector(k,n)
Range = n-k+2;
for i = k+1:n+1
    t(i) = i-k;
end
t(n+2:n+k+1) = Range;
t = t./Range;   % normalize
end

% First Order Base Function
function N1 = FirstOrderBSplineFunction(k,t,x,SampleSize)
m1 = length(t); m2 = length(x);
N1 = zeros(m1-1,m2);
k1 = k;
k2 = m1-k;
for i = k1:k2
    j1 = (i-1)*SampleSize+1;
    j2 = j1+SampleSize-1;
    N1(i,j1:j2) = 1;
end
end

% General Order Base Function
function N = GeneralOrderBsplineFunction(k,t,x,SampleSize,N1)
m1 = length(t); 
N = N1;
for l = 2:k
    d = l;
    k1 = k-(d-1);
    k2 = m1-k;

    i = k1;
    j1 = (i-1)*SampleSize+1;
    j2 = j1+d*SampleSize-1;
    N(i,j1:j2) = (t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
    for i = k1+1:k2-1
        j1 = (i-1)*SampleSize+1;
        j2 = j1+d*SampleSize-1;
        N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
    end
    i = k2;
    j1 = (i-1)*SampleSize+1;
    j2 = j1+d*SampleSize-1;
    N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2);
end
end