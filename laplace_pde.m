
%KANE NORMAN

%PDE'S AND FINITE DIFFERENCE METHODS
%PLEASE READ WRITE UP!


  

%#######PART 1 JACOBI ITERARTION#########

 p=1;%INDEX FOR SUBPLOT GRAPHS
 err=zeros(1,3); %STORES ERRORS FOR OUR 3 ITERATIVE METHODS
 entry=1;%COUNTER USED TO INDEX TABLE 2 AND ERROR
 Table1=zeros(3,4);%PREALLOCATE TABLE FOR PART 1
 Table2=zeros(3,4);%PREALLOCATE TABLE FOR PART 2

 
for h=[1/4,1/8,1/16]
 
    n=(1/h)-1; %SET SIZE N
    tol=10^(-8);%STOPPING CRITERIA
    ratio=1;%INSTANTIATE RATIO
    k=0;%COUNTER FOR NUM ITERATIONS
    
   
    
    w_true=zeros(n+2,n+2); 
%CONSTRUCT A MATRIX OF OUR TRUE SOLUTION
    for i=1:n+2
        for j=1:n+2
            w_true(i,j)=sin(pi*(j-1)*h)*exp(pi*(i-1)*h); 
        end
    end

    
  %CONSTRUCT A MATRIX FOR OUR INITAL GUESS W 
  w_old=zeros(n+2,n+2); %ASSIGN THE INTERIOR VALUES THE VALUE ZERO
    for i=1:n+2 %ASSIGN THE BOUNDARY CONDITIONS
        for j=1:n+2
    w_old(end,j)=exp(pi)*sin(pi*(j-1)*h);
    w_old(1,j)=sin(pi*(j-1)*h);
    
    w_old(i,end)=0;
    w_old(i,1)=0;
   
        end
    end
    

    w=zeros(size(w_old));
    for i=1:n+2  %ASSIGN THE BOUNDARY CONDITIONS TO W
        for j=1:n+2
            w(end,j)=exp(pi)*sin(pi*(j-1)*h);
            w(1,j)=sin(pi*(j-1)*h);
            w(i,end)=0;
            w(i,1)=0;
        end
    end
    
       
    %INITIALIZE f AND residual
    f=zeros(size(w));
    res=zeros(size(f));


 while ratio>tol
   k=k+1; %UPDATE COUNTER/ITERATION NUMBER
   
    %CALCULATE AN APPROXIMATION FOR f
    for i=2:n+1
        for j=2:n+1
           res(i,j)=(4*w_old(i,j)-w_old(i+1,j)-w_old(i-1,j)-w_old(i,j+1)-w_old(i,j-1))/(h^2);
        end
    end
    
    %COMPUTE RESIDUAL
    r(k)=max(max(abs(f-res)));
    %COMPUTE RATIO
    ratio=r(k)/r(1);
  
    %USE JACOBIS METHOD TO UPDATE OUR INITIAL GUESS
    for i=2:n+1
        for j=2:n+1
            w(i,j)=((w_old(i+1,j)+w_old(i-1,j)+w_old(i,j-1)+w_old(i,j+1))+f(i,j)*h^2)/4; 
        end
    end
    
   %UPDATE Wk to Wk+1 
   w_old=w;
   
 end
 
%PLOT CONTOUR MAP AND MESH GRID
 x=0:h:1;
 y=0:h:1;
 [x,y]=meshgrid(x,y);
  %CONTOUR PLOT
 subplot(3,2,p);
 axis square
 contour(x,y,w)
 xlabel("x position")
 ylabel("y position")
 str=sprintf('h=%s',strtrim(rats(h)));
 title(str);
 %MESH PLOT
 subplot(3,2,p+1)
 mesh(x,y,w)
 xlabel("x position")
 ylabel("y position")
 zlabel("temperature")
 str=sprintf('h=%s',strtrim(rats(h)));
 title(str);
 %UPDATE P BY TWO SO ALL 6 FIGURES FIT IN THE SUBPLOT
 p=p+2;
 
 err(entry)=max(max(abs(w-w_true)));
 
 
 %CREATING TABLE FOR ANALYSIS OF JACOBI'S METHOD
 Table1(entry,1)=h;
 Table1(entry,2)=err(entry);
 
 if entry>1
 Table1(entry,3)=err(entry-1)/err(entry);
 Table1(entry,4)=log2(err(entry-1)/err(entry));
 end
 
 
Table2(entry,2)=k;
entry=entry+1;
end
 Table1=array2table(Table1,'VariableNames',{'h','error','ratio','order'}) %DISPLAY RESULTS FROM PART 1 IN A TABLE

 

%******PART 2 GS ITERATION*******

    %%WE WILL USE THE SAME METHODS AS IN PART 1. EXCEPT WE WILL USE THE G-S ITERATIVE METHOD

entry=1;
for h=[1/4,1/8,1/16]
 
    n=(1/h)-1; %SET n
    tol=10^(-8);%STOPPING CRITERIA 
    ratio=1;%INITALIZE RATIO
    k=0;%INITALIZE COUNTING VARIABLE 
    
   
    
    w_true=zeros(n+2,n+2); 
%CONSTRUCT A MATRIX OF OUR TRUE SOLUTION
    for i=1:n+2
        for j=1:n+2
            w_true(i,j)=sin(pi*(j-1)*h)*exp(pi*(i-1)*h); 
        end
    end

    
 
    w=zeros(size(n+2,n+2));
    for i=1:n+2  %ASSIGN THE BOUNDARY CONDITIONS TO W
        for j=1:n+2
            w(end,j)=exp(pi)*sin(pi*(j-1)*h);
            w(1,j)=sin(pi*(j-1)*h);
            w(i,end)=0;
            w(i,1)=0;
        end
    end
    
       
    f=zeros(size(w));
    res=zeros(size(f));


 while ratio>tol
   k=k+1; %UPDATE COUNTER/NUM ITERATIONS
   
    %APPROXIMATE F TO CALCULATE RESIDUAL
    for i=2:n+1
        for j=2:n+1
           res(i,j)=(4*w(i,j)-w(i+1,j)-w(i-1,j)-w(i,j+1)-w(i,j-1))/(h^2);
        end
    end
    
    %COMPUTING RATIO
    r(k)=max(max(abs(f-res)));
    ratio=r(k)/r(1);
  
    %USE GS METHOD TO UPDATE OUR INITAL GUESS
    for i=2:n+1
        for j=2:n+1
            w(i,j)=((w(i+1,j)+w(i-1,j)+w(i,j-1)+w(i,j+1))+f(i,j)*h^2)/4; 
        end
    end
    
 
 end
 Table2(entry,3)=k;
 entry=entry+1;
end



% %**********PART 3 SOR*********************
entry=1;

for h=[1/4,1/8,1/16]
 
    n=(1/h)-1; %SET SIZE OF N
    tol=10^(-8);%STOPPING CRITERIA
    ratio=1;%INITIALIZE RATIO
    k=0;%COUNTER VARIABLE
    
    omega=(2)/(1+sqrt(1-(cos(pi*h))^2));
   
    
    w_true=zeros(n+2,n+2); 
%CONSTRUCT THE MATRIX OF OUR TRUE SOLUTION phi(x,y)=sin(pi*x)*exp(pi*y)
    for i=1:n+2
        for j=1:n+2
            w_true(i,j)=sin(pi*(j-1)*h)*exp(pi*(i-1)*h); 
        end
    end

    
 
    w=zeros(size(n+2,n+2));
    for i=1:n+2  %ASSIGN THE BOUNDARY CONDITIONS TO W
        for j=1:n+2
            w(end,j)=exp(pi)*sin(pi*(j-1)*h);
            w(1,j)=sin(pi*(j-1)*h);
            w(i,end)=0;
            w(i,1)=0;
        end
    end
    
    
    f=zeros(size(w));
    res=zeros(size(f));


 while ratio>tol
   k=k+1; %COUNTER FOR NUMBER OF ITERATIONS
   
    %CALCULATE AN APPROXIMATION TO F
    for i=2:n+1
        for j=2:n+1
           res(i,j)=(4*w(i,j)-w(i+1,j)-w(i-1,j)-w(i,j+1)-w(i,j-1))/(h^2);
        end
    end
    
    %COMPUTE RESIDUAL
    r(k)=max(max(abs(f-res)));
    %COMPUTE RATIO
    ratio=r(k)/r(1);
  
    %USE SOR TO UPDATE OUR INITIAL GUESS
    for i=2:n+1
        for j=2:n+1
            w(i,j)=w(i,j)+(omega/4)*(f(i,j)*h^2 +w(i+1,j)+w(i-1,j)+w(i,j+1)+w(i,j-1)-4*w(i,j));
        end
    end
    

 end
 
 Table2(entry,4)=k;
 entry=entry+1;

end

%TABLE FOR COMPARING THE 3 ITERATIVE METHODS
Table2(1,1)=1/4;
Table2(2,1)=1/8;
Table2(3,1)=1/64;
Table2 = array2table(Table2,'VariableNames',{'h','Jacobi','GS','SOR'}) 