clc;clear 
sum = zeros(1,1536)
    K=1;
    N1=512;    
    M1=512;      
    Phi=PartFourierMtx(M1,N1)/sqrt(N1); 
    a = Gaussian();
    x = [a(1:512)].'
    s=Phi*x;                                        
    m=2*K;                                            
    Psi=dwtmtx(512,'coif2',1)./sqrt(N1);                        
    T=Phi*Psi';                                      
    hat_y=zeros(1,N1);                                                    
    Aug_t1=[];                                         
    r_n=s;                                            
    for times=1:1;                                    
        for col=1:N1;                                  
            product1(col)=abs(T(:,col)'*r_n);          
        end
        [val,pos]=max(product1);                      
        Aug_t1=[Aug_t1,T(:,pos)];                       
    
        T(:,pos)=zeros(M1,1);                          
        aug_y=(Aug_t1'*Aug_t1)^(-1)*Aug_t1'*s;           
        r_n=s-Aug_t1*aug_y;                            
        pos_array(times)=pos;                         
    end
    hat_y(pos_array)=aug_y;
    hat_x=abs(real(Psi'*hat_y.'));                         
    c =  Gaussian()
    d = c(1:1536)
    x_r = d.';
    A = hat_x;
    N2=512;    
    M2=512;      
    Phi=PartFourierMtx(M2,N2)/sqrt(N2); 
    a =  Gaussian()
    x2 = [a(513:1024)].';
    s2=Phi*x2;                                        
    m=2*K;                                            
    Psi=dwtmtx(512,'coif2',1)./sqrt(N2);                        
    T2=Phi*Psi';                                      
    hat_y=zeros(1,N2);                                                    
    Aug_t2=[];                                         
    r_n=s2;                                            
    for times=1:1;                                    
        for col=1:N2;                                  
            product2(col)=abs(T2(:,col)'*r_n);          
        end
        [val,pos]=max(product2);                      
        Aug_t2=[Aug_t2,T2(:,pos)];                       
    
        T2(:,pos)=zeros(M2,1);                          
        aug_y=(Aug_t2'*Aug_t2)^(-1)*Aug_t2'*s2;           
        r_n=s2-Aug_t2*aug_y;                            
        pos_array(times)=pos;                         
    end
    hat_y(pos_array)=aug_y;
    hat_x2=abs(real(Psi'*hat_y.'))
    N3=512;    
    M3=512;      
    Phi=PartFourierMtx(M3,N3)/sqrt(N3); 
    a =  Gaussian()
    x3 = [a(1025:1536)].';
    s3=Phi*x3;                                                                                 
    Psi=dwtmtx(512,'coif2',1)./sqrt(N3);                        
    T3=Phi*Psi';                                      
    hat_y=zeros(1,N3);                                                    
    Aug_t3=[];                                         
    r_n=s3;                                            
    for times=1:1;                                    
        for col=1:N3;                                  
            product3(col)=abs(T3(:,col)'*r_n);          
        end
        [val,pos]=max(product3);                      
        Aug_t3=[Aug_t3,T3(:,pos)];                       
    
        T3(:,pos)=zeros(M3,1);                          
        aug_y=(Aug_t3'*Aug_t3)^(-1)*Aug_t3'*s3;           
        r_n=s3-Aug_t3*aug_y;                            
        pos_array(times)=pos;                           
    end
    hat_y(pos_array)=aug_y;

    hat_x3=abs(real(Psi'*hat_y.'))
    A2=hat_x3;
    figure(1);
    B = [hat_x.',hat_x2.',hat_x3.']
hold on;
C = 250:1:1785
plot(C,B)
plot(C,x_r,'r');
%legend('Recovery','Original')
xlswrite('C:\Users\fangzheng\Desktop\1.xlsx',B)
