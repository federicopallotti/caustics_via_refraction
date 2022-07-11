A = imread('test3.png')/256; %input image, with values scaled down to be from 0 to 1
%input image resolution scaling
numrows = 10;
numcols = 10;
C = imresize(A,[numrows numcols]);
imshow(C)
C=double(C);
mapping_u = zeros(numrows,numcols);
mapping_v = zeros(numrows,numcols);

%light term
L = sum(C,'all')/(numrows*numcols);

% initialize the mapping
for i=1:numrows
    for j=1:numcols
           mapping_u(i,j) = i;
           mapping_v(i,j) = j;
    end
end
mapping = cat(3,double(mapping_u),double(mapping_v));
%delaunay triangulation of the grid, for the sake of visualization
DT = delaunay(mapping_u, mapping_v);

% main execution
delta = 0.1; %derivative step
alpha = 0.01; %learning rate
d=10; %distance between the 2 screens
m = main(mapping,delta,alpha,numrows,numcols,L,C,DT);
norms = compute_norms(m,numrows,numcols,d);
zs = zeros(numrows,numcols);
[x,y] = meshgrid(1:numcols,1:numrows);
quiver3(x,y,zs, norms(:,:,1),norms(:,:,2),norms(:,:,3));
axis equal

% write a txt file with X,Y,Z  U,V,W (coordinates and components of each
% normal)
[x,y] = meshgrid(1:numcols,1:numrows);
first = x(:);
second = y(:);
third = zeros(length(first),1);

N_x = norms(:,:,1);
N_y = norms(:,:,2);
N_z = norms(:,:,3);

T = table(first, second,third, N_x(:), N_y(:), N_z(:), 'VariableNames', { 'X', 'Y','Z','U','V','W'} );
writetable(T, 'MyFile.txt');
compute_angles(norms);
%visualize results as an image 
vis = visualize(numcols,numrows,m);
imshow(vis);

%function to visualize the results
function img_res = visualize(numcols,numrows,map)
    tmp = zeros(numcols,numrows);
    f_x = map(:,:,1);
    f_y = map(:,:,2);
    for i=1:numrows-1
        for j=1:numcols-1
            area = (f_x(i+1,j) -f_x(i,j)) * (f_y(i,j+1) -f_y(i,j));
            tmp(i,j) = (area/(numcols*numrows))*255;        
        end
    end
    img_res = tmp;
end

% function for angles computation from components (used for surface
% construction on blender)
function angles = compute_angles(norms)
    [theta,rho,z] =  cart2pol(norms(:,:,1),norms(:,:,2),norms(:,:,3));
    theta = rad2deg(theta);
    rho = rad2deg(rho);
    angles = cat(3,theta,rho);
    
end

% function that returns the matrix of normals
function normals = compute_norms(map,numrows,numcols,d)
    normals = zeros(numrows,numcols,3);
    mu = 1.52; %refraction ratio
    in = [0; 0; -1]; %incoming light ray
    
    for i=1:numrows
        for j=1:numcols
               t = [map(i,j,1) - i; map(i,j,2) - (numrows - j); -d];
               t = t/norm(t); %normalized 
               n = (t-mu*in)/norm(t-mu*in);
               normals(i,j,1) = n(1);
               normals(i,j,2) = n(2);
               normals(i,j,3) = n(3);
        end
    end
    %   interpolation
    [x,y] = meshgrid(1:numcols,1:numrows);
    
    first = map(x(:),y(:),1);
    second = map(x(:),y(:),2);
    third = normals(x(:),y(:),1);
    
    N_x_F = scatteredInterpolant(first(:),second(:),third(:));
    
    third = normals(x(:),y(:),2);
    
    N_y_F = scatteredInterpolant(first(:),second(:),third(:));
    
    third = normals(x(:),y(:),3);
    
    N_z_F = scatteredInterpolant(first(:),second(:),third(:));
    
    N = zeros(numcols,numrows,3);
    N(:,:,1) = N_x_F(x,y);
    N(:,:,2) = N_y_F(x,y);
    N(:,:,3) = N_z_F(x,y);
  
    normals = N;
    
    
end


% main function 
function m = main(map,delta,alpha,numrows,numcols,L,C, DT)
    previous_error = 100000000;
    while (alpha > 0.00001 && delta > 0.00001)
          G = gradient(map,delta,numrows,numcols,L,C);
          new_map = zeros(numrows,numcols,2);
          
%         loop for checking if actually increment position of a point
            f_x = map(:,:,1);
            f_y = map(:,:,2);
            for i=1:numrows
                for j=1:numcols
                    if i==1 && j==1
                        a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                        b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
                        
                        det1 = det([a;b]);
                        
                        if det1>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end

                    elseif  i == numcols && j==1
                       
                        b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];

                        c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];

                        det3 = det([b;c]);

                        if   det3>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end

                    elseif i == numcols && j == numrows
                        c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
                        d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];
                        det2 = det([c;d]);
                        
                        if det2>0 && det3>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                    
                    elseif i ==1 && j == numrows
                        a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                        d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];
                        det4 = det([d;a]);

                        if det4>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                        
                    elseif i ==1
                        a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                        b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
                        d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];

                        det1 = det([a;b]);
                        det4 = det([d;a]);
                        
                        if det1>0 && det4>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                        
                    elseif j ==1
                        a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                        b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
                        c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
                        
                        det1 = det([a;b]);
                        det3 = det([b;c]);
                        
                        if det1>0 && det3>0 
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                        
                    elseif i == numcols
                        
                        b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
                        c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
                        d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];

                        det2 = det([c;d]);
                        det3 = det([b;c]);
                        
                        if  det2>0 && det3>0 
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                    elseif j ==numrows
                        a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                        c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
                        d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];
 
                        det2 = det([c;d]);
                        det4 = det([d;a]);
                        
                        if det2>0 && det4>0
                            new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                        else
                            new_map(i,j,:) = map(i,j,:); 
                        end
                        
                    else
                        
                    a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
                    b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
                    c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
                    d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];

                    det1 = det([a;b]);
                    det2 = det([c;d]);
                    det3 = det([b;c]);
                    det4 = det([d;a]);
                    
                    if det1>0 && det2>0 && det3>0 && det4>0
                        new_map(i,j,:) = bound(map(i,j,:) - alpha*G(i,j,:),1,numcols);
                    else
                        new_map(i,j,:) = map(i,j,:); 
                    end
                    end
                      
                end
            end
          
        
          new_map(1,:,1) = map(1,:,1); 
          new_map(numcols,:,1) = map(numcols,:,1); 
          new_map(:,1,2) = map(:,1,2); 
          new_map(:,numcols,2) = map(:,numcols,2);
         
        plot(new_map(:,:,1),new_map(:,:,2),'*');
        triplot(DT,new_map(:,:,1),new_map(:,:,2));
        %uncomment to visualize the error function
%       E = abs(error_matrix(new_map,numcols,numrows,L,C))/200;
        %imshow(E);
        %E(1,1)
        drawnow;
        current_error = error(new_map,numcols,numrows,L,C);
        
        if previous_error > current_error 
            map = new_map;
            previous_error = current_error;
        else
            alpha = alpha/2;
        end
        
        fprintf("Error %f,  stepsize, %f \n", current_error, alpha);
    end
    
    m = map;
    
end

% function to compute the gradient
function G = gradient(map,delta,numrows,numcols,L,C)
% gradient matrix
    G_x = zeros(numrows,numcols);
    G_y = zeros(numrows,numcols);
    
    for i=1:numrows
        for j=1:numcols
%             once for the y
              map_new_y = update_map(map,delta,'y',i,j);
              map = double(map);
              
              G_y(i,j) = (error(map_new_y, numcols,numrows,L,C) - error(map,numcols,numrows,L,C))/delta;
%             once for the x
              map_new_x = update_map(map,delta,'x',i,j);
              G_x(i,j) = (error(map_new_x, numcols,numrows,L,C) - error(map,numcols,numrows,L,C))/delta;
               
        end
    end
    
    G = cat(3,double(G_x),double(G_y));
    
end


% function to compute the error for one of the N component of the gradient
function e = error(map,numcols,numrows,L,C)
    Jacobian = zeros(numrows,numcols);
    for i=1:numrows
        for j=1:numcols
           Jacobian(i,j) = 0.25*jacobian_el(i,j,map,numcols,numrows);

        end
    end
    % compute the D
    D = L*Jacobian - C;
    e = sqrt(sum(D.^2,'all'));
    
end

%function to compute error matrix
function D = error_matrix(map,numcols,numrows,L,C)
    Jacobian = zeros(numrows,numcols);
    for i=1:numrows
        for j=1:numcols
           Jacobian(i,j) = jacobian_el(i,j,map,numcols,numrows);

        end
    end
    
    % compute the D
    D = L*Jacobian - C;
   
end

% function to compute the el at pos i,j of the jacobian matrix
function determinant = jacobian_el(i,j,map,numcols,numrows)
    f_x = map(:,:,1);
    f_y = map(:,:,2);

    if i==1
        i = i+1;
    end
    
    if j==1
        j = j+1;
    end
    
    if i == numcols
        i = i-1;
    end
        
    if j == numrows
        j = j-1;
    end
    
    a = [(f_x(i+1,j)-f_x(i,j)),(f_y(i+1,j)-f_y(i,j))];
    b = [(f_x(i,j+1)-f_x(i,j)),(f_y(i,j+1)-f_y(i,j))];
    c = [f_x(i-1,j)-f_x(i,j), f_y(i-1,j)-f_y(i,j)];
    d = [f_x(i,j-1)-f_x(i,j),f_y(i,j-1)-f_y(i,j)];
    
    det1 = det([a;b]);
    det2 = det([c;d]);
    det3 = det([b;c]);
    det4 = det([d;a]);
    
    determinant = det1*0.25 + det2*0.25 + det3*0.25 + det4*0.25;
    
end

% func to compute the mapping with a specific displacement
function m = update_map(map,delta,component,i,j)
    f_x = map(:,:,1);
    f_y = map(:,:,2);
    if component=='x'
        f_x(i,j) = f_x(i,j) + delta;
    end
    
    if component == 'y'
        f_y(i,j) = f_y(i,j) + delta;
    end
    
    m = cat(3,double(f_x),double(f_y));
        
end


function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
   y=min(max(x,bl),bu);
end






    
       
           
            